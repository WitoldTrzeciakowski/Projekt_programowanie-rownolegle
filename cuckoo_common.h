#ifndef CUCKOO_COMMON_H
#define CUCKOO_COMMON_H

#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <random>
#include <thread>

using namespace std;

// ===== RNG =====
thread_local static mt19937 rng(random_device{}() + hash<thread::id>{}(this_thread::get_id()));
thread_local static uniform_real_distribution<double> dist(0.0, 1.0);

// ===== helpers =====
inline double randDouble(const double a, const double b) {
    return a + (b - a) * ((double)rand() / RAND_MAX);
}

inline double randDoubleParallel(const double a, const double b) {
    return a + (b - a) * dist(rng);
}

inline double levyFlight() {
    double u = randDouble(0, 1);
    double v = randDouble(0, 1);
    return u / pow(fabs(v), 1.0 / 1.5);
}

inline double levyFlightParallel() {
    double u = dist(rng);
    double v = dist(rng);
    return u / pow(fabs(v), 1.0 / 1.5);
}

// ===== interfaces =====
vector<double> cuckooSearch(
    int, int, int, double, double, double,
    const function<double(const vector<double>&)>
&);

vector<double> cuckooSearchParallel(
    int, int, int, double, double, double,
    function<double(const vector<double>&)>, int = 0
);

vector<double> cuckooSearchProcess(
    int, int, int, double, double, double,
    const function<double(const vector<double>&)>&, int = 0
);

#endif
