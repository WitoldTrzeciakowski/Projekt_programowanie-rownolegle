#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <algorithm>

using namespace std;

// ===================== helpers =====================


double randDouble(double a, double b) {
    return a + (b - a) * ((double)rand() / RAND_MAX);
}


double levyFlight() {
    double u = randDouble(0, 1);
    double v = randDouble(0, 1);
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
