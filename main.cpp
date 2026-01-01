#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include "helper_functions.cpp"
double sphereFunction(const std::vector<double>& x) {
    double sum = 0.0;
    for (double v : x)
        sum += v * v;
    return sum;
}
int main() {
    srand(time(NULL));

    vector<double> best = cuckooSearch(
        25,     // liczba gniazd
        2,      // wymiar
        500,    // iteracje
        0.25,   // pa
        -10.0,  // dolna granica
        10.0,   // górna granica
        sphereFunction // funkcja celu
    );

    cout << "Najlepsze rozwiazanie:\n";
    for (double x : best)
        cout << x << " ";
    cout << endl;

    cout << "Wartosc funkcji: " << sphereFunction(best) << endl;
}
