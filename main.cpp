#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <string>
#include <cmath>
#include <chrono>

#include "helper_functions.cpp"

using namespace std;
using namespace chrono;

// ===================== TEST FUNCTIONS =====================

double sphereFunction(const vector<double>& x) {
    double sum = 0.0;
    for (double v : x)
        sum += v * v;
    return sum;
}

double rosenbrockFunction(const vector<double>& x) {
    double sum = 0.0;
    for (size_t i = 0; i < x.size() - 1; i++)
        sum += 100 * pow(x[i+1] - x[i]*x[i], 2) + pow(x[i] - 1, 2);
    return sum;
}

double rastriginFunction(const vector<double>& x) {
    double sum = 10.0 * x.size();
    for (double v : x)
        sum += v*v - 10.0 * cos(2 * M_PI * v);
    return sum;
}

double ackleyFunction(const vector<double>& x) {
    double sum1 = 0.0, sum2 = 0.0;
    for (double v : x) {
        sum1 += v * v;
        sum2 += cos(2 * M_PI * v);
    }
    return -20.0 * exp(-0.2 * sqrt(sum1 / x.size()))
           - exp(sum2 / x.size()) + 20 + M_E;
}

double schwefelFunction(const vector<double>& x) {
    double sum = 0.0;
    for (double v : x)
        sum += v * sin(sqrt(abs(v)));
    return 418.9829 * x.size() - sum;
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

// ===================== RUN + TIME =====================

vector<double> runTimed(
    const string& name,
    function<vector<double>()> algo,
    function<double(const vector<double>&)> fitness
) {
    auto start = high_resolution_clock::now();
    vector<double> result = algo();
    auto end = high_resolution_clock::now();

    double time_ms = duration<double, milli>(end - start).count();

    cout << "--------------------------------------\n";
    cout << "Wersja: " << name << "\n";
    cout << "Czas:   " << time_ms << " ms\n";
    cout << "f(x):   " << fitness(result) << "\n";
    cout << "x*:     ";
    for (double v : result) cout << v << " ";
    cout << "\n";

    return result;
}

// ===================== MAIN =====================

int main(int argc, char* argv[]) {
    srand(time(NULL));

    // ----- default parameters -----
    int numNests   = 5000;
    int iterations = 30;
    int dimension  = 500;
    string func    = "sphere";
    int numProc    = 100;

    // ----- CLI override -----
    if (argc > 1) numNests   = fromString(argv[1], numNests);
    if (argc > 2) iterations = fromString(argv[2], iterations);
    if (argc > 3) dimension  = fromString(argv[3], dimension);
    if (argc > 4) func       = argv[4];
    if (argc > 5) numProc    = fromString(argv[5], numProc);

    // ----- bounds & fitness -----
    double LB, UB;
    function<double(const vector<double>&)> fitness;

    if (func == "rosen") {
        fitness = rosenbrockFunction;
        LB = -5.0; UB = 10.0;
    }
    else if (func == "rastrigin") {
        fitness = rastriginFunction;
        LB = -5.12; UB = 5.12;
    }
    else if (func == "ackley") {
        fitness = ackleyFunction;
        LB = -32.0; UB = 32.0;
    }
    else if (func == "schwefel") {
        fitness = schwefelFunction;
        LB = -500.0; UB = 500.0;
    }
    else {
        fitness = sphereFunction;
        LB = -10.0; UB = 10.0;
    }

    const double pa = 0.25;

    cout << "======================================\n";
    cout << "Cuckoo Search – porownanie czasow\n";
    cout << "Funkcja: " << func << "\n";
    cout << "Wymiar:  " << dimension << "\n";
    cout << "Nests:   " << numNests << "\n";
    cout << "Iter:    " << iterations << "\n";
    cout << "======================================\n";

    // ===================== RUN ALL =====================

    runTimed("SEQUENTIAL", [&]() {
        return cuckooSearch(
            numNests, dimension, iterations,
            pa, LB, UB, fitness
        );
    }, fitness);

    runTimed("PARALLEL (threads)", [&]() {
        return cuckooSearchParallel(
            numNests, dimension, iterations,
            pa, LB, UB, fitness
        );
    }, fitness);

    runTimed("PROCESS", [&]() {
        return cuckooSearchProcess(
            numNests, dimension, iterations,
            pa, LB, UB, fitness, numProc
        );
    }, fitness);

    return 0;
}
