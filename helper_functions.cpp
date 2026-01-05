#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <thread>
#include <random>

using namespace std;

// ===================== helpers =====================

// Thread-local random number generator
thread_local mt19937 rng(random_device{}() + hash<thread::id>{}(this_thread::get_id()));
thread_local uniform_real_distribution<double> dist(0.0, 1.0);

double randDouble(double a, double b) {
    return a + (b - a) * ((double)rand() / RAND_MAX);
}

double levyFlight() {
    double u = randDouble(0, 1);
    double v = randDouble(0, 1);
    return u / pow(fabs(v), 1.0 / 1.5);
}

double randDoubleParrarel(double a, double b) {
    return a + (b - a) * dist(rng);
}

double levyFlightParrarel() {
    double u = dist(rng);
    double v = dist(rng);
    return u / pow(fabs(v), 1.0 / 1.5);
}

vector<double> cuckooSearch(int N, int DIM, int MAX_ITER,
                            double pa, double LB, double UB,
                            function<double(const vector<double>&)> fitness)
{
    vector<vector<double>> nests(N, vector<double>(DIM));
    vector<double> fitnessVal(N);

    // Inicjalizacja
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < DIM; j++)
            nests[i][j] = randDouble(LB, UB);
        fitnessVal[i] = fitness(nests[i]);
    }

    for (int iter = 0; iter < MAX_ITER; iter++) {

        // --- Lévy flights ---
        for (int i = 0; i < N; i++) {
            vector<double> candidate = nests[i];

            for (int d = 0; d < DIM; d++) {
                candidate[d] += levyFlight();
                candidate[d] = max(LB, min(UB, candidate[d]));
            }

            double f_new = fitness(candidate);
            if (f_new < fitnessVal[i]) {
                nests[i] = candidate;
                fitnessVal[i] = f_new;
            }
        }

        // --- Porzucanie gniazd ---
        for (int i = 0; i < N; i++) {
            if ((double)rand() / RAND_MAX < pa) {
                for (int d = 0; d < DIM; d++)
                    nests[i][d] = randDouble(LB, UB);
                fitnessVal[i] = fitness(nests[i]);
            }
        }
    }

    // Wybór najlepszego rozwiązania
    int best = min_element(fitnessVal.begin(), fitnessVal.end()) - fitnessVal.begin();
    return nests[best];
}

// ===================== PARALLEL VERSION =====================

void levyFlightWorker(vector<vector<double>>& nests, 
                      vector<double>& fitnessVal,
                      int start, int end, int DIM, 
                      double LB, double UB,
                      function<double(const vector<double>&)> fitness)
{
    for (int i = start; i < end; i++) {
        vector<double> candidate = nests[i];

        for (int d = 0; d < DIM; d++) {
            candidate[d] += levyFlightParrarel();
            candidate[d] = max(LB, min(UB, candidate[d]));
        }

        double f_new = fitness(candidate);
        if (f_new < fitnessVal[i]) {
            nests[i] = candidate;
            fitnessVal[i] = f_new;
        }
    }
}

void abandonNestsWorker(vector<vector<double>>& nests,
                        vector<double>& fitnessVal,
                        int start, int end, int DIM,
                        double pa, double LB, double UB,
                        function<double(const vector<double>&)> fitness)
{
    for (int i = start; i < end; i++) {
        if (randDoubleParrarel(0, 1) < pa) {
            for (int d = 0; d < DIM; d++)
                nests[i][d] = randDoubleParrarel(LB, UB);
            fitnessVal[i] = fitness(nests[i]);
        }
    }
}

vector<double> cuckooSearchParallel(int N, int DIM, int MAX_ITER,
                                     double pa, double LB, double UB,
                                     function<double(const vector<double>&)> fitness,
                                     int num_threads = 0)
{
    if (num_threads <= 0)
        num_threads = thread::hardware_concurrency();
    
    vector<vector<double>> nests(N, vector<double>(DIM));
    vector<double> fitnessVal(N);

    // Inicjalizacja (parallel)
    vector<thread> threads;
    int chunk = N / num_threads;
    
    for (int t = 0; t < num_threads; t++) {
        int start = t * chunk;
        int end = (t == num_threads - 1) ? N : (t + 1) * chunk;
        
        threads.emplace_back([&, start, end]() {
            for (int i = start; i < end; i++) {
                for (int j = 0; j < DIM; j++)
                    nests[i][j] = randDoubleParrarel(LB, UB);
                fitnessVal[i] = fitness(nests[i]);
            }
        });
    }
    
    for (auto& t : threads)
        t.join();
    threads.clear();

    // Main loop
    for (int iter = 0; iter < MAX_ITER; iter++) {

        // --- Lévy flights (parallel) ---
        for (int t = 0; t < num_threads; t++) {
            int start = t * chunk;
            int end = (t == num_threads - 1) ? N : (t + 1) * chunk;
            
            threads.emplace_back(levyFlightWorker, ref(nests), ref(fitnessVal),
                               start, end, DIM, LB, UB, fitness);
        }
        
        for (auto& t : threads)
            t.join();
        threads.clear();

        // --- Porzucanie gniazd (parallel) ---
        for (int t = 0; t < num_threads; t++) {
            int start = t * chunk;
            int end = (t == num_threads - 1) ? N : (t + 1) * chunk;
            
            threads.emplace_back(abandonNestsWorker, ref(nests), ref(fitnessVal),
                               start, end, DIM, pa, LB, UB, fitness);
        }
        
        for (auto& t : threads)
            t.join();
        threads.clear();
    }

    // Wybór najlepszego rozwiązania
    int best = min_element(fitnessVal.begin(), fitnessVal.end()) - fitnessVal.begin();
    return nests[best];
}
// ===================== PARALLEL VERSION – PROCESSES =====================
// Linux / Unix only

#include <sys/mman.h>
#include <sys/wait.h>
#include <unistd.h>

// ---------- shared memory helpers ----------

double* allocSharedDouble(size_t n) {
    return (double*) mmap(
        nullptr,
        n * sizeof(double),
        PROT_READ | PROT_WRITE,
        MAP_SHARED | MAP_ANONYMOUS,
        -1,
        0
    );
}

inline double& NEST(double* nests, int i, int d, int DIM) {
    return nests[i * DIM + d];
}

// ---------- worker: Lévy flights ----------

void levyFlightProcess(double* nests, double* fitnessVal,
                       int start, int end, int DIM,
                       double LB, double UB,
                       function<double(const vector<double>&)> fitness)
{
    vector<double> candidate(DIM);

    for (int i = start; i < end; i++) {
        for (int d = 0; d < DIM; d++) {
            candidate[d] = NEST(nests, i, d, DIM) + levyFlight();
            candidate[d] = max(LB, min(UB, candidate[d]));
        }

        double f_new = fitness(candidate);
        if (f_new < fitnessVal[i]) {
            for (int d = 0; d < DIM; d++)
                NEST(nests, i, d, DIM) = candidate[d];
            fitnessVal[i] = f_new;
        }
    }
}

// ---------- worker: abandon nests ----------

void abandonNestsProcess(double* nests, double* fitnessVal,
                         int start, int end, int DIM,
                         double pa, double LB, double UB,
                         function<double(const vector<double>&)> fitness)
{
    vector<double> tmp(DIM);

    for (int i = start; i < end; i++) {
        if (randDouble(0, 1) < pa) {
            for (int d = 0; d < DIM; d++) {
                tmp[d] = randDouble(LB, UB);
                NEST(nests, i, d, DIM) = tmp[d];
            }
            fitnessVal[i] = fitness(tmp);
        }
    }
}

// ---------- MAIN: cuckoo search (process-based) ----------

vector<double> cuckooSearchProcess(int N, int DIM, int MAX_ITER,
                                   double pa, double LB, double UB,
                                   function<double(const vector<double>&)> fitness,
                                   int num_proc)
{
    if (num_proc <= 0)
        num_proc = sysconf(_SC_NPROCESSORS_ONLN);

    double* nests = allocSharedDouble(N * DIM);
    double* fitnessVal = allocSharedDouble(N);

    // ----- initialization -----
    for (int i = 0; i < N; i++) {
        vector<double> tmp(DIM);
        for (int d = 0; d < DIM; d++) {
            tmp[d] = randDouble(LB, UB);
            NEST(nests, i, d, DIM) = tmp[d];
        }
        fitnessVal[i] = fitness(tmp);
    }

    int chunk = N / num_proc;

    for (int iter = 0; iter < MAX_ITER; iter++) {

        // ----- Lévy flights -----
        for (int p = 0; p < num_proc; p++) {
            if (fork() == 0) {
                int start = p * chunk;
                int end = (p == num_proc - 1) ? N : (p + 1) * chunk;
                levyFlightProcess(nests, fitnessVal,
                                  start, end, DIM, LB, UB, fitness);
                _exit(0);
            }
        }
        while (wait(nullptr) > 0);

        // ----- abandon nests -----
        for (int p = 0; p < num_proc; p++) {
            if (fork() == 0) {
                int start = p * chunk;
                int end = (p == num_proc - 1) ? N : (p + 1) * chunk;
                abandonNestsProcess(nests, fitnessVal,
                                    start, end, DIM, pa, LB, UB, fitness);
                _exit(0);
            }
        }
        while (wait(nullptr) > 0);
    }

    // ----- best solution -----
    int best = 0;
    for (int i = 1; i < N; i++)
        if (fitnessVal[i] < fitnessVal[best])
            best = i;

    vector<double> result(DIM);
    for (int d = 0; d < DIM; d++)
        result[d] = NEST(nests, best, d, DIM);

    munmap(nests, N * DIM * sizeof(double));
    munmap(fitnessVal, N * sizeof(double));

    return result;
}
