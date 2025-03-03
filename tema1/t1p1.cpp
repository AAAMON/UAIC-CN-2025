#include <iostream>


void update_u(float& u)
{
	u = u/10.0;
}


int main ()
{
	std::cout << "Hello World\n";

	int m = 0;
	float prev_u = -1;
	float u = 1;

	// MAIN LOOP
	while (!(1 + u == 1))
	{
		++m;
		prev_u = u;
		update_u(u);
		std::cout << "u is: " << u << '\n';
	}

	std::cout << "Precizia masina e: " << prev_u << '\n';
	// resulting u is 1.4013e-45

	return 0;
}
