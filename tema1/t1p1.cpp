// Task: search for the smallest number, before rounding errors make the sum insignificant/pointless.

#include <iostream>

void update_u(double& u) // u is of form: 10^(-m)
{
	u = u/10.0;
}

int main ()
{
	int m = 0;
	double prev_u = -1;
	double u = 1;

	// MAIN LOOP
	while (!(1 + u == 1))
	{
		++m;
		prev_u = u;
		update_u(u);
		std::cout << "u este: " << u << '\n';
	}

	std::cout << "Precizia masina e: " << prev_u << '\n';
	// Andy: resulting u is 1.4013e-45 (originally)
	// Ama: u = 1e-15

	return 0;
}