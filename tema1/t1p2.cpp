// Task: Check associativity based on results from previous exercise.
// Task: Find a set of numbers where the multiplication is not associative.

#include <iostream>

int main ()
{
	// we took this value from the previous exercise (see file t1p1.cpp)
	double u = 1e-15; // ((1 + u) != 1) AND ((1 + u/10) == 1) are true

	double x = 1.0;
	double y = u/10; // 1e-16
	double z = u/10; // 1e-16

	std::cout << "Relatia (x + y) + z != x + (y + z) este evaluata ca: " << ((x + y) + z != x + (y + z)) << '\n';
    // (x + y) + z = 1
    // x + (y + z) > 1

    double a = 1e15; // the opposite of b
    double b = 1e-15; // same value as u, because of c
    // a * b = 1 (aprox.)
    double c = 1 + u; // should be != 1

    std::cout << "Relatia (a * b) * c != a * (b * c) este evaluata ca: " << ((a * b) * c != a * (b * c)) << '\n';
    // (a * b) * c = c
    // a * (b * c) = a * (u + u^2) = a * b + a * u^2 = 1 + a * u^2 > 1 + u = c

	return 0;
}
