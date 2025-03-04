// Task: Check associativity based on results from previous exercise.
// Task: Find a set of numbers where the multiplication is not associative.

#include <iostream>

int main ()
{
	// we took this value from the previous exercise (see file t1p1.cpp)
	float u = 1e-08;

	float x = 1.0;
	float y = u/10.0;
	float z = u/10.0;

	std::cout << "Relatia (x + y) + z != x + (y + z) este evaluata ca: " << ((x + y) + z != x + (y + z)) << '\n';

	return 0;
}
