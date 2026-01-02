#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <string>
#include "helper_functions.cpp"
double sphereFunction(const std::vector<double>& x) {
    double sum = 0.0;
    for (double v : x)
        sum += v * v;
    return sum;
}

// Helper to convert string to type safely
template <typename T>
T fromString(const char* str, T defaultValue) {
    std::istringstream ss(str);
    T value;
    if (ss >> value)
        return value;
    return defaultValue;
}

int main(int argc, char* argv[]) {
    srand(time(NULL));

    // Default parameters
    int numNests = 300000;
    int iterations = 1000;
    string version = "seq"; // default sequential

    // Override from command line if provided
    if (argc > 1) numNests   = fromString(argv[1], numNests);
    if (argc > 2) iterations = fromString(argv[2], iterations);
    if (argc > 3) version    = argv[3];

    // Hardcoded parameters
    const int dimension = 2;
    const double pa = 0.25;
    const double lowerBound = -10.0;
    const double upperBound = 10.0;

    vector<double> best;

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
    } else {
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

    cout << "Najlepsze rozwiazanie (" << version << "):\n";
    for (double x : best)
        cout << x << " ";
    cout << endl;

    cout << "Wartosc funkcji: " << sphereFunction(best) << endl;
}
