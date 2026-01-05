#include "cuckoo_common.h"

vector<double> cuckooSearch(int N, int DIM, int MAX_ITER,
                            double pa, double LB, double UB,
                            function<double(const vector<double>&)> fitness)
{
    vector<vector<double>> nests(N, vector<double>(DIM));
    vector<double> fitnessVal(N);

    for (int i = 0; i < N; i++) {
        for (int d = 0; d < DIM; d++)
            nests[i][d] = randDouble(LB, UB);
        fitnessVal[i] = fitness(nests[i]);
    }

    for (int iter = 0; iter < MAX_ITER; iter++) {

        for (int i = 0; i < N; i++) {
            vector<double> candidate = nests[i];
            for (int d = 0; d < DIM; d++) {
                candidate[d] += levyFlight();
                candidate[d] = max(LB, min(UB, candidate[d]));
            }
            double f = fitness(candidate);
            if (f < fitnessVal[i]) {
                nests[i] = candidate;
                fitnessVal[i] = f;
            }
        }

        for (int i = 0; i < N; i++) {
            if (randDouble(0, 1) < pa) {
                for (int d = 0; d < DIM; d++)
                    nests[i][d] = randDouble(LB, UB);
                fitnessVal[i] = fitness(nests[i]);
            }
        }
    }

    int best = min_element(fitnessVal.begin(), fitnessVal.end()) - fitnessVal.begin();
    return nests[best];
}
