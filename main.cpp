#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <string>
#include <cmath>
#include <chrono>
#include <functional>
#include <fstream>

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

struct TimedResult {
    vector<double> solution;
    double fx;
    double time_ms;
};

TimedResult runTimed(
    const string& name,
    function<vector<double>()> algo,
    function<double(const vector<double>&)> fitness
) {
    auto start = high_resolution_clock::now();
    vector<double> result = algo();
    auto end = high_resolution_clock::now();

    double time_ms = duration<double, milli>(end - start).count();
    double fx = fitness(result);

    cout << "--------------------------------------\n";
    cout << "Wersja: " << name << "\n";
    cout << "Czas:   " << time_ms << " ms\n";
    cout << "f(x):   " << fx << "\n";
    cout << "x*:     ";
    for (double v : result) cout << v << " ";
    cout << "\n";

    return {result, fx, time_ms};
}

// ===================== BENCHMARK STRUCT =====================

struct Benchmark {
    string name;
    function<double(const vector<double>&)> fitness;
    double LB;
    double UB;
};

// ===================== RUN ONE BENCHMARK =====================

void runBenchmark(
    const Benchmark& bm,
    int numNests,
    int iterations,
    int dimension,
    int numProc,
    ofstream* csv = nullptr
) {
    const double pa = 0.25;

    cout << "\n======================================\n";
    cout << "Funkcja: " << bm.name << "\n";
    cout << "Wymiar:  " << dimension << "\n";
    cout << "Nests:   " << numNests << "\n";
    cout << "Iter:    " << iterations << "\n";
    cout << "======================================\n";

    vector<pair<string, TimedResult>> results;

    // results.push_back({"SEQUENTIAL", runTimed("SEQUENTIAL", [&]() {
    //     return cuckooSearch(numNests, dimension, iterations, pa, bm.LB, bm.UB, bm.fitness);
    // }, bm.fitness)});

    results.push_back({"THREADS", runTimed("THREADS", [&]() {
        return cuckooSearchParallel(numNests, dimension, iterations, pa, bm.LB, bm.UB, bm.fitness, numProc);
    }, bm.fitness)});

    results.push_back({"PROCESS", runTimed("PROCESS", [&]() {
        return cuckooSearchProcess(numNests, dimension, iterations, pa, bm.LB, bm.UB, bm.fitness, numProc);
    }, bm.fitness)});

    // Write to CSV if provided
    if (csv && csv->is_open()) {
        for (auto& r : results) {
            string numProcStr = (r.first == "PROCESS") ? to_string(numProc) : "";
            *csv << bm.name << ","        // Function
                 << r.first << ","        // Method
                 << numNests << ","       // Nests
                 << iterations << ","     // Iterations
                 << dimension << ","      // Dimension
                 << numProcStr << ","     // NumProc (only for PROCESS)
                 << r.second.time_ms      // Time_ms
                 << "\n";
        }
    }
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
    string csvFile = "";

    // ----- CLI override -----
    if (argc > 1) numNests   = fromString(argv[1], numNests);
    if (argc > 2) iterations = fromString(argv[2], iterations);
    if (argc > 3) dimension  = fromString(argv[3], dimension);
    if (argc > 4) func       = argv[4];
    if (argc > 5) numProc    = fromString(argv[5], numProc);
    if (argc > 6) csvFile    = argv[6];

    // ----- CSV file -----
    ofstream csv;
    bool writeHeader = true;

    if (!csvFile.empty()) {
        // Check if file exists
        ifstream check(csvFile);
        if (check.good()) {
            writeHeader = false; // file exists → do not write header
        }
        check.close();

        csv.open(csvFile, ios::app); // append mode
        if (!csv.is_open()) {
            cerr << "Nie mozna otworzyc pliku CSV: " << csvFile << "\n";
            return 1;
        }
        if (writeHeader)
            csv << "Function,Method,Nests,Iterations,Dimension,NumProc,Time_ms\n";
    }

    // ----- benchmarks -----
    vector<Benchmark> benchmarks = {
        {"sphere",    sphereFunction,    -10.0,   10.0},
        {"rosen",     rosenbrockFunction, -5.0,    10.0},
        {"rastrigin", rastriginFunction,  -5.12,   5.12},
        {"ackley",    ackleyFunction,     -32.0,   32.0},
        {"schwefel",  schwefelFunction,   -500.0,  500.0}
    };

    // ----- run -----
    if (func == "all") {
        for (const auto& bm : benchmarks) {
            runBenchmark(bm, numNests, iterations, dimension, numProc, csvFile.empty() ? nullptr : &csv);
        }
    } else {
        bool found = false;
        for (const auto& bm : benchmarks) {
            if (bm.name == func) {
                runBenchmark(bm, numNests, iterations, dimension, numProc, csvFile.empty() ? nullptr : &csv);
                found = true;
                break;
            }
        }
        if (!found) {
            cerr << "Unknown function: " << func << "\n";
            return 1;
        }
    }

    if (csv.is_open()) csv.close();

    return 0;
}
