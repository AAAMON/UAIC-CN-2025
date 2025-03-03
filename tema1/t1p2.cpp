#include <iostream>

int main ()
{
	std::cout << "Hello World\n";

	// we took this value from the previous exercise. See file t1p1.cpp
	float u = 1e-07;

	float x = 1.0;
	float y = u/10.0;
	float z = u/10.0;

	std::cout << "Relatia (x+y)+z != x+(y+z) este evaluata ca: " << ((x+y)+z != x + (y+z));

	return 0;
}
