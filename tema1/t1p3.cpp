// Task: Polinome hierarchy based on approximations of the sine function

#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
// pt sort
#include <algorithm>
#include <chrono>

std::ofstream file("t1p3_output.txt");
std::ifstream infile("t1p3_output.txt");
// constants:
double c1 = 0.16666666666666666666666666666667;
double c2 = 0.00833333333333333333333333333333;
double c3 = 1.984126984126984126984126984127e-4;
double c4 = 2.7557319223985890652557319223986e-6;
double c5 = 2.5052108385441718775052108385442e-8;
double c6 = 1.6059043836821614599392377170155e-10;

void calculeaza_polinoame_si_scrie_top3(std::ofstream &file, double victim, std::vector<int> times)
{
    double correct_val = std::sin(victim);

    double help1 = victim * victim;
    // help1 e x^2

    using namespace std::chrono;
    auto start = high_resolution_clock::now();
    double P1 = victim * (1 + help1 * (-c1 + c2 * help1));
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start).count();
    times[0] = duration;

    double help2 = help1 * help1 * help1 * victim;
    // help2 e x^7

    start = high_resolution_clock::now();
    double P2 = P1 - c3 * help2;
    end = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(end-start).count();
    times[1] = duration;

    double help3 = -c3 * help2 + c4 * help2 * help1;
    // help3 e -c3*x^7 + c4*x^9

    start = high_resolution_clock::now();
    double P3 = P1 + help3;
    end = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(end-start).count();
    times[2] = duration;

    start = high_resolution_clock::now();
    double P4 = victim * (1 + help1 * (-0.166 + help1 * 0.00833)) + help3;
    end = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(end-start).count();
    times[3] = duration;

    start = high_resolution_clock::now();
    double P5 = victim * (1 + help1 * (-0.1666 + help1 * 0.008333)) + help3;
    end = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(end-start).count();
    times[4] = duration;

    start = high_resolution_clock::now();
    double P6 = victim * (1 + help1 * (-0.16666 + help1 * 0.0083333)) + help3;
    end = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(end-start).count();
    times[5] = duration;

    double help4 = help2 * help1 * help1;
    // help4 e x^11

    start = high_resolution_clock::now();
    double P7 = P3 - c5 * help4;
    end = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(end-start).count();
    times[6] = duration;

    start = high_resolution_clock::now();
    double P8 = P7 + c6 * help4 * help1;
    end = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(end-start).count();
    times[7] = duration;

    // calculam direct si erorile
    std::vector<std::pair<double, std::string>> polinoamee =
        {
            {std::abs(P1 - correct_val), "1"},
            {std::abs(P2 - correct_val), "2"},
            {std::abs(P3 - correct_val), "3"},
            {std::abs(P4 - correct_val), "4"},
            {std::abs(P5 - correct_val), "5"},
            {std::abs(P6 - correct_val), "6"},
            {std::abs(P7 - correct_val), "7"},
            {std::abs(P8 - correct_val), "8"},
        };
    std::sort(polinoamee.begin(), polinoamee.end());

    file << victim << '-' << polinoamee[0].second << ' ' << polinoamee[1].second << ' ' << polinoamee[2].second
         << ' ' << polinoamee[3].second << ' ' << polinoamee[4].second << '\n';
}

int main()
{

    // scriem top 3 polinoame pt fiecare nr intr-un fisier bc it's like 30k...
    if (!file)
    {
        std::cerr << "Nu se deschide fisierul bruh\n";
        return 1;
    }

    const int NUMBERS_COUNT = 10000;
    // Pe romaneste, daca dai print la numar o sa fie aprox. intre -1.57 si 1.57
    const double LOWER_BOUND = -M_PI / 2;
    const double UPPER_BOUND = M_PI / 2;

    // facem setup la functia random
    std::random_device rd;
    std::mt19937 gen(rd()); // Alg. Mersenne Twister
    std::uniform_real_distribution<double> dist(LOWER_BOUND, UPPER_BOUND);

    double current_victim;

    using namespace std::chrono;
    auto start = high_resolution_clock::now();

    std::vector<int> times = {0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<std::pair<int, std::string>> times_median =
    {
        {0, "1"}, 
        {0, "2"}, 
        {0, "3"}, 
        {0, "4"}, 
        {0, "5"}, 
        {0, "6"}, 
        {0, "7"}, 
        {0, "8"}
    };    

    // MAIN LOOP
    for (int i = 1; i <= NUMBERS_COUNT; ++i)
    {
        // generam random folosind distributia dist() cu generatorul gen
        current_victim = dist(gen);
        std::cout << "#" << i << ": Number: " << current_victim << '\n';
        calculeaza_polinoame_si_scrie_top3(file, current_victim, times);

        for(int j = 0; j < 8; j++)
        {
            std::cout<<"Polinomul #"<<j<<": "<<times[j]<<"(ms)\n";
            times_median[j].first += times[j];
        }
        std::cout<<"----------\n";
    }
    for(int i = 0; i < 8; i++)
    {
        times_median[i].first /= NUMBERS_COUNT;
    }
    std::sort(times_median.begin(), times_median.end());

    if (!infile)
    {
        std::cerr << "Nu se deschide fisierul bruh\n";
        return 1;
    }

    std::vector<std::pair<int, std::string>> scores = 
        {
            {0, "1"}, 
            {0, "2"}, 
            {0, "3"}, 
            {0, "4"}, 
            {0, "5"}, 
            {0, "6"}, 
            {0, "7"}, 
            {0, "8"}
        };
    std::string line;
    while(std::getline(infile, line))
    {
        if(line.empty())
        {
            continue;
        }
        // sine values can be negative, but all have a floating point anyway, so we start looking for the '-' from position 1 (not 0)
        int pos = line.find('-', 1);
        std::istringstream iss(line.substr(pos + 1)); // takes everything after the line
        for(int i = 0; i < 5; i++) // we ranked the 5 best in the output file
        {
            int poly;
            iss >> poly;
            scores[poly - 1].first += 5 - i; // first place gets the highest score
        }
    }
    std::sort(scores.begin(), scores.end());
    std::reverse(scores.begin(), scores.end());
    std::cout << "Ierarhia polinoamelor: \n";
    for(int i = 0; i < 8; i++)
    {
        std::cout << "#" << i + 1 << ": " << scores[i].second << "(scor: " << scores[i].first << ") \n";
    }
    std::cout<<"Timpii de lucru: \n";
    for(int i = 0; i < 8; i++)
    {
        std::cout << "#" << i + 1 << ": " << times_median[i].first << "(ms) \n";
    }

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start).count();

    std::cout << "Toata nebunia pentru " << NUMBERS_COUNT << " iteratii a durat: " << duration << "ms\n";

    return 0;
}
