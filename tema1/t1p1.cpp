#include <iostream>


void update_u(float& u)
{
	u = u/10.0;
}


int main ()
{
	std::cout << "Hello World\n";

	int m = 0;
	float u = 1;

	// MAIN LOOP
	while (!(1 + u = 1))
	{
		++m;
		update_u(u);
		std::cout << "u is: " << u << '\n';
	}

	return 0;
}
