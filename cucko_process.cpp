#include "cuckoo_common.h"
#include <sys/mman.h>
#include <sys/wait.h>
#include <unistd.h>
#include <semaphore.h>

static double* allocShared(const size_t n) {
    return (double*) mmap(nullptr, n * sizeof(double),
        PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
}

// Współdzielone semafory
struct SharedData {
    sem_t start_sem;
    sem_t done_sem;
    int iter;
    bool stop;
};

vector<double> cuckooSearchProcess(const int N, const int DIM, const int MAX_ITER,
                                   const double pa, const double LB, const double UB,
                                   const function<double(const vector<double>&)>& fitness,
                                   int num_proc)
{
    if (num_proc <= 0)
        num_proc = sysconf(_SC_NPROCESSORS_ONLN);

    double* nests = allocShared(N * DIM);
    double* fit = allocShared(N);
    
    // Współdzielona struktura kontrolna
    SharedData* shared = (SharedData*)mmap(nullptr, sizeof(SharedData),
        PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
    
    sem_init(&shared->start_sem, 1, 0);
    sem_init(&shared->done_sem, 1, 0);
    shared->iter = 0;
    shared->stop = false;

    // Inicjalizacja populacji
    for (int i = 0; i < N; i++) {
        vector<double> tmp(DIM);
        for (int d = 0; d < DIM; d++) {
            tmp[d] = randDouble(LB, UB);
            nests[i * DIM + d] = tmp[d];
        }
        fit[i] = fitness(tmp);
    }

    // Tworzenie procesów potomnych
    for (int p = 0; p < num_proc; p++) {
        if (fork() == 0) {
            // Proces potomny
            int chunk = N / num_proc;
            int s = p * chunk;
            int e = (p == num_proc - 1) ? N : (p + 1) * chunk;
            
            while (true) {
                // Czekaj na start iteracji
                sem_wait(&shared->start_sem);
                
                if (shared->stop) {
                    _exit(0);
                }
                
                // Wykonaj swoją część pracy
                int current_iter = shared->iter;
                
                for (int i = s; i < e; i++) {
                    // Generowanie nowego rozwiązania
                    vector<double> cand(DIM);
                    for (int d = 0; d < DIM; d++) {
                        cand[d] = nests[i * DIM + d] + levyFlight();
                        if (cand[d] < LB) cand[d] = LB;
                        if (cand[d] > UB) cand[d] = UB;
                    }
                    
                    double f = fitness(cand);
                    
                    // Atomowe porównanie (jeśli potrzebne)
                    if (f < fit[i]) {
                        for (int d = 0; d < DIM; d++) {
                            nests[i * DIM + d] = cand[d];
                        }
                        __sync_synchronize(); // Bariera pamięci
                        fit[i] = f;
                    }
                }
                
                // Powiadom o zakończeniu
                sem_post(&shared->done_sem);
            }
        }
    }

    // Główna pętla w rodzicu
    for (int iter = 0; iter < MAX_ITER; iter++) {
        shared->iter = iter;
        
        // Rozpocznij iterację - obudź wszystkie procesy
        for (int p = 0; p < num_proc; p++) {
            sem_post(&shared->start_sem);
        }
        
        // Czekaj na zakończenie wszystkich procesów
        for (int p = 0; p < num_proc; p++) {
            sem_wait(&shared->done_sem);
        }
        
        // Abandon worst nests (sekwencyjnie, bo mało pracy)
        for (int i = 0; i < N; i++) {
            if ((double)rand() / RAND_MAX < pa) {
                for (int d = 0; d < DIM; d++) {
                    nests[i * DIM + d] = randDouble(LB, UB);
                }
                fit[i] = fitness(vector<double>(nests + i * DIM, 
                                              nests + i * DIM + DIM));
            }
        }
    }
    
    // Zakończ procesy potomne
    shared->stop = true;
    for (int p = 0; p < num_proc; p++) {
        sem_post(&shared->start_sem);
    }
    
    // Czekaj na procesy potomne
    while (wait(nullptr) > 0);

    // Znajdź najlepsze rozwiązanie
    int best = 0;
    for (int i = 1; i < N; i++) {
        if (fit[i] < fit[best]) best = i;
    }

    vector<double> res(DIM);
    for (int d = 0; d < DIM; d++) {
        res[d] = nests[best * DIM + d];
    }

    // Zwolnij zasoby
    munmap(nests, N * DIM * sizeof(double));
    munmap(fit, N * sizeof(double));
    sem_destroy(&shared->start_sem);
    sem_destroy(&shared->done_sem);
    munmap(shared, sizeof(SharedData));

    return res;
}