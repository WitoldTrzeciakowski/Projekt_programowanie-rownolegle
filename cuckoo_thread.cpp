#include "cuckoo_common.h"

static void levyWorker(vector<vector<double>>& nests,
                       vector<double>& fitnessVal,
                       int s, int e, int DIM,
                       double LB, double UB,
                       function<double(const vector<double>&)> fitness)
{
    for (int i = s; i < e; i++) {
        vector<double> cand = nests[i];
        for (int d = 0; d < DIM; d++) {
            cand[d] += levyFlightParallel();
            cand[d] = max(LB, min(UB, cand[d]));
        }
        double f = fitness(cand);
        if (f < fitnessVal[i]) {
            nests[i] = cand;
            fitnessVal[i] = f;
        }
    }
}

static void abandonWorker(vector<vector<double>>& nests,
                          vector<double>& fitnessVal,
                          int s, int e, int DIM,
                          double pa, double LB, double UB,
                          function<double(const vector<double>&)> fitness)
{
    for (int i = s; i < e; i++) {
        if (randDoubleParallel(0, 1) < pa) {
            for (int d = 0; d < DIM; d++)
                nests[i][d] = randDoubleParallel(LB, UB);
            fitnessVal[i] = fitness(nests[i]);
        }
    }
}

vector<double> cuckooSearchParallel(int N, int DIM, int MAX_ITER,
                                    double pa, double LB, double UB,
                                    function<double(const vector<double>&)> fitness,
                                    int num_threads)
{
    if (num_threads <= 0)
        num_threads = thread::hardware_concurrency();

    vector<vector<double>> nests(N, vector<double>(DIM));
    vector<double> fitnessVal(N);

    for (int i = 0; i < N; i++) {
        for (int d = 0; d < DIM; d++)
            nests[i][d] = randDoubleParallel(LB, UB);
        fitnessVal[i] = fitness(nests[i]);
    }

    int chunk = N / num_threads;

    for (int iter = 0; iter < MAX_ITER; iter++) {

        vector<thread> threads;
        for (int t = 0; t < num_threads; t++) {
            int s = t * chunk;
            int e = (t == num_threads - 1) ? N : (t + 1) * chunk;
            threads.emplace_back(levyWorker,
                                 ref(nests), ref(fitnessVal),
                                 s, e, DIM, LB, UB, fitness);
        }
        for (auto& th : threads) th.join();

        threads.clear();
        for (int t = 0; t < num_threads; t++) {
            int s = t * chunk;
            int e = (t == num_threads - 1) ? N : (t + 1) * chunk;
            threads.emplace_back(abandonWorker,
                                 ref(nests), ref(fitnessVal),
                                 s, e, DIM, pa, LB, UB, fitness);
        }
        for (auto& th : threads) th.join();
    }

    int best = min_element(fitnessVal.begin(), fitnessVal.end()) - fitnessVal.begin();
    return nests[best];
}
