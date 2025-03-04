// Task: Polinome hierarchy based on approximations of the sine function

#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <fstream>
// pt sort
#include <algorithm>

std::ofstream file("t1p3_output.txt");
// constants:
double c1 = 0.16666666666666666666666666666667;
double c2 = 0.00833333333333333333333333333333;
double c3 = 1.984126984126984126984126984127e-4;
double c4 = 2.7557319223985890652557319223986e-6;
double c5 = 2.5052108385441718775052108385442e-8;
double c6 = 1.6059043836821614599392377170155e-10;

double P1(double x)
{
    double y = x * x;
    return x * (1 + y * (- c1 + c2 * y));
}

double P2(double x)
{
    double y = x * x;
    return x * (1 + y * (- c1 + y * (c2 - c3 * y)));
}

double P3(double x)
{
    double y = x * x;
    return x * (1 + y * (- c1 + y * (c2 + y * (- c3  + c4 * y))));
}

double P4(double x)
{
    double y = x * x;
    return x * (1 + y * (-0.166 + y * (0.00833 + y * (- c3  + c4 * y))));
}

double P5(double x)
{
    double y = x * x;
    return x * (1 + y * (-0.1666 + y * (0.008333 + y * (- c3  + c4 * y))));
}

double P6(double x)
{
    double y = x * x;
    return x * (1 + y * (-0.16666 + y * (0.0083333 + y * (- c3  + c4 * y))));
}

double P7(double x)
{
    double y = x * x;
    return x * (1 + y * (- c1 + y * (c2 + y * (- c3  + y * (c4 - c5 * y)))));
}

double P8(double x)
{
    double y = x * x;
    return x * (1 + y * (- c1 + y * (c2 + y * (- c3  + y * (c4 + y * (- c5 + c6 * y))))));
}

void calculeaza_polinoame_si_scrie_top3(std::ofstream& file, double victim)
{
    double correct_val = std::sin(victim);

    // calculam direct si erorile
    std::vector<std::pair<double, std::string>> polinoamee =
    {
	{std::abs(P1(victim) - correct_val), "1"},
	{std::abs(P2(victim) - correct_val), "2"},
	{std::abs(P3(victim) - correct_val), "3"},
	{std::abs(P4(victim) - correct_val), "4"},
	{std::abs(P5(victim) - correct_val), "5"},
	{std::abs(P6(victim) - correct_val), "6"},
	{std::abs(P7(victim) - correct_val), "7"},
	{std::abs(P8(victim) - correct_val), "8"},
    };
    std::sort(polinoamee.begin(), polinoamee.end());
    
    file << victim << '-' << polinoamee[0].second << ' ' << polinoamee[1].second << ' ' << polinoamee[2].second
    << ' ' << polinoamee[3].second << ' ' << polinoamee[4].second << '\n';
}

int main ()
{
    // scriem top 3 polinoame pt fiecare nr intr-un fisier bc it's like 30k...
    
    if (!file)
    {
        std::cerr << "Nu se deschide fisierul bruh\n";
	    return 1;
    }

    const int NUMBERS_COUNT = 100;
    // Pe romaneste, daca dai print la numar o sa fie aprox. intre -1.57 si 1.57
    const double LOWER_BOUND = -M_PI / 2; 
    const double UPPER_BOUND = M_PI / 2;

    // facem setup la functia random
    std::random_device rd;
    std::mt19937 gen(rd()); // Alg. Mersenne Twister
    std::uniform_real_distribution<double> dist(LOWER_BOUND, UPPER_BOUND);

    double current_victim;

    // MAIN LOOP
    for (int i = 1; i <= NUMBERS_COUNT; ++i)
    {
	    // generam random folosind distributia dist() cu generatorul gen
	    current_victim = dist(gen);
        std::cout << "#" << i << ": Number: " << current_victim << '\n';
        calculeaza_polinoame_si_scrie_top3(file, current_victim);
    }

    // citim datele din fisier si facem o ierarhie

    return 0;
}
