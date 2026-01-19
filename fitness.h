#ifndef FITNESS_H
#define FITNESS_H

#include <functional>
#include <vector>
#include "benchmark.h"
#include "shared_types.h"

inline std::function<double(const std::vector<double>&)> getFitness(int id) {
    switch (id) {
    case FITNESS_SPHERE:     return sphereFunction;
    case FITNESS_ROSENBROCK: return rosenbrockFunction;
    case FITNESS_RASTRIGIN:  return rastriginFunction;
    case FITNESS_ACKLEY:     return ackleyFunction;
    case FITNESS_SCHWEFEL:   return schwefelFunction;
    default:                 return sphereFunction;
    }
}

#endif //FITNESS_H