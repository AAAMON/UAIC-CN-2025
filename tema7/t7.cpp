#include <iostream>
#include <vector>
#include <cmath> // pentru fabs()
#include <fstream>
#include <cstdlib> // pentru rand()
#include <random>
using namespace std;

// Functie pentru calculul intervalului [-R, R]
double calculeaza_R(const vector<double> &coeficienti)
{
    double a0 = coeficienti[0];
    double A = 0.0;

    // Calculam A = max(|a1|, ..., |an|)
    for (int i = 1; i < coeficienti.size(); ++i)
    {
        A = max(A, fabs(coeficienti[i]));
    }

    // Calculam R folosind formula
    double R = (fabs(a0) + A) / fabs(a0);
    return R;
}

// Evaluare Horner si derivate
double horner_eval(const vector<double> &coef, double x)
{
    double result = coef[0];
    for (int i = 1; i < coef.size(); ++i)
    {
        result = result * x + coef[i];
    }
    return result;
}

// Derivata de ordin 1
vector<double> der1(const vector<double> &coef)
{
    vector<double> deriv;
    int grad = coef.size() - 1; // Ignora termenul liber (constantele derivate dau mereu 0)
    for (int i = 0; i < grad; ++i)
    {
        deriv.push_back(coef[i] * (grad - i)); // Explicatie: P(x) = 1 * x^3 − 6 * x^2 + 11 * x − 6 => P'(x) = 3 * x^2 − 12 * x + 11
                                               // 1 * (3 - 0) * x^2 - 6 * (3 - 1) * x^1 + 11 * (3 - 2) * x, ultimul coeficient fiind sarit automat
    }
    return deriv;
}

// Derivata de ordin 2
vector<double> der2(const vector<double> &coef)
{
    return der1(der1(coef));
}

// Verifica daca o radacina este distincta
bool distinct(double r, const vector<double> &rad, const double epsilon)
{
    for (double x : rad)
    {
        if (fabs(x - r) < epsilon)
        {
            return false;
        }
    }
    return true;
}

// Metoda lui Halley
bool halley(const vector<double> &coef, double x0, double &radacina, const double epsilon, const int kmax)
{
    vector<double> d1 = der1(coef);
    vector<double> d2 = der2(coef);
    double x = x0; // x0 se alege aleator (vezi functia main)
    int k = 0;

    while (k < kmax) // In loc de do-while (k = 0)
    {
        double Px = horner_eval(coef, x); // P(xk)
        double P1x = horner_eval(d1, x);  // P'(xk)
        double P2x = horner_eval(d2, x);  // P''(xk)

        double A = 2 * P1x * P1x - Px * P2x; // numaratorul pt ak
        if (fabs(A) < epsilon)
        {
            cout << "Impartire la 0" << endl;
            return false; // A prea mic, risc de impartire la 0
        }

        double delta = 2 * Px * P1x / A;
        x = x - delta;
        ++k;

        if (fabs(delta) < epsilon)
        {
            radacina = x;
            return true;
        }

        if (fabs(delta) > 1e8)
        {
            cout << "Divergenta" << endl;
            return false; // Divergenta
        }
    }

    cout << "Depasit numarul maxim de pasi" << endl;
    return false; // Am depasit numarul de pasi
}

int main()
{
    const double epsilon = 1e-6; // Precizia
    const int kmax = 1000;       // Numar maxim de iteratii
    const int x0_cnt = 100;      // Cate puncte x0 sa testam

    // Exemplu de polinom: x^3 - 6 * x^2 + 11 * x - 6
    // vector<double> coef = {1.0, -6.0, 11.0, -6.0}; // a0, a1, ..., an

    // Exemplu de polinom: (1 / 8) * (8 * x^4 - 38 * x^3 + 49 * x^2 - 22 * x + 3)
    vector<double> coef = {8.0, -38.0, 49.0, -22.0, 3.0}; // a0, a1, ..., an
    for (double &c : coef)
    {
        c /= 8.0; // normalizare
    }

    // Exemplu de polinom: x^4 - 6 * x^3 + 13 * x^2 - 12 * x + 4
    // vector<double> coef = {1.0, -6.0, 13.0, -12.0, 4.0}; // a0, a1, ..., an
    // Rezultatul asteptat: 1, 2
    // Rezultatul real: Impartire la 0. Singura solutie: micsoram epsilon. Problema? Va fi mult prea mic epsilon, si toate radacinile vor fi considerate 'distincte'.

    double R = calculeaza_R(coef);
    cout << "Toate radacinile reale se afla in intervalul [-" << R << ", " << R << "]" << endl;

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(-R, R); // x0 in [-R, R]

    vector<double> radacini;

    for (int i = 0; i < x0_cnt; ++i)
    {
        double x0 = dist(gen); // Alegem aleator x0 in [-R, R], pentru a fi in vecinatatea lui x*
        double r = 0.0;

        if (halley(coef, x0, r, epsilon, kmax))
        {
            if (distinct(r, radacini, epsilon))
            {
                radacini.push_back(r);
                cout << "Radacina gasita: " << r << "\n";
            }
        }
    }

    // Scriem in fisier
    ofstream fout("radacini.txt");
    for (double r : radacini)
    {
        fout << r << "\n";
    }
    fout.close();

    cout << "Radacinile distincte au fost salvate in fisierul 'radacini.txt'.\n";
    return 0;
}
