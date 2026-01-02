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