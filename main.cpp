#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <string>

#include "helper_functions.cpp"

using namespace std;

// ===================== TEST FUNCTION =====================

double sphereFunction(const vector<double>& x) {
    double sum = 0.0;
    for (double v : x)
        sum += v * v;
    return sum;
}

// ===================== STRING → TYPE =====================

template <typename T>
T fromString(const char* str, T defaultValue) {
    istringstream ss(str);
    T value;
    if (ss >> value)
        return value;
    return defaultValue;
}

// ===================== MAIN =====================

int main(int argc, char* argv[]) {
    srand(time(NULL));

    // ----- default parameters -----
    int numNests   = 300;
    int iterations = 100;
    string version = "seq";   // seq | par | proc
    int numProc    = 0;       // only for proc

    // ----- CLI override -----
    if (argc > 1) numNests   = fromString(argv[1], numNests);
    if (argc > 2) iterations = fromString(argv[2], iterations);
    if (argc > 3) version    = argv[3];
    if (argc > 4) numProc    = fromString(argv[4], numProc);

    // ----- problem parameters -----
    const int dimension = 2;
    const double pa = 0.25;
    const double lowerBound = -10.0;
    const double upperBound = 10.0;

    vector<double> best;

    // ===================== SELECT VERSION =====================

    if (version == "par") {
        best = cuckooSearchParallel(
            numNests,
            dimension,
            iterations,
            pa,
            lowerBound,
            upperBound,
            sphereFunction
        );
    }
    else if (version == "proc") {
        best = cuckooSearchProcess(
            numNests,
            dimension,
            iterations,
            pa,
            lowerBound,
            upperBound,
            sphereFunction,
            numProc   // 0 = auto
        );
    }
    else { // seq
        best = cuckooSearch(
            numNests,
            dimension,
            iterations,
            pa,
            lowerBound,
            upperBound,
            sphereFunction
        );
    }

    // ===================== OUTPUT =====================

    cout << "Najlepsze rozwiazanie (" << version << "):\n";
    for (double x : best)
        cout << x << " ";
    cout << endl;

    cout << "Wartosc funkcji: " << sphereFunction(best) << endl;

    return 0;
}
