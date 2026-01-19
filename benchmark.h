#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <vector>
#include <string>
#include <cmath>
#include <functional>

// ===================== TEST FUNCTIONS =====================

inline double sphereFunction(const std::vector<double>& x) {
    double sum = 0.0;
    for (const double v : x)
        sum += v * v;
    return sum;
}

inline double rosenbrockFunction(const std::vector<double>& x) {
    double sum = 0.0;
    for (size_t i = 0; i < x.size() - 1; i++)
        sum += 100 * pow(x[i+1] - x[i]*x[i], 2) + pow(x[i] - 1, 2);
    return sum;
}

inline double rastriginFunction(const std::vector<double>& x) {
    const auto size = static_cast<double>(x.size());
    double sum = 10.0 * size;
    for (const double v : x)
        sum += v*v - 10.0 * cos(2 * M_PI * v);
    return sum;
}

inline double ackleyFunction(const std::vector<double>& x) {
    double sum1 = 0.0, sum2 = 0.0;
    for (const double v : x) {
        sum1 += v * v;
        sum2 += cos(2 * M_PI * v);
    }
    const auto size = static_cast<double>(x.size());
    return -20.0 * exp(-0.2 * sqrt(sum1 / size))
           - exp(sum2 / size) + 20 + M_E;
}

inline double schwefelFunction(const std::vector<double>& x) {
    const auto size = static_cast<double>(x.size());
    double sum = 0.0;
    for (const double v : x)
        sum += v * sin(sqrt(abs(v)));
    return 418.9829 * size - sum;
}


// ===================== BENCHMARK STRUCT =====================

struct Benchmark {
    std::string name;
    int fitnessId;
    double LB;
    double UB;
};


#endif //BENCHMARK_H