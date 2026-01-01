#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include "helper_functions.cpp"
int main() {
    srand(time(NULL));

    vector<double> best = cuckooSearch(
        25,     // liczba gniazd
        2,      // wymiar
        500,    // iteracje
        0.25,   // pa
        -10.0,  // dolna granica
        10.0    // górna granica
    );

    cout << "Najlepsze rozwiazanie:\n";
    for (double x : best)
        cout << x << " ";
    cout << endl;

    cout << "Wartosc funkcji: " << fitness(best) << endl;
}
