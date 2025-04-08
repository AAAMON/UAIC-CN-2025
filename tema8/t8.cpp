#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <functional>

// ---------------- Config ----------------
const double epsilon = 1e-4;
const int kmax = 30000;
const double h = 1e-5;
const int n_dim = 2; // numarul de variabile
const bool backtracking = true;
// ----------------------------------------

// Functia de test
double F(const std::vector<double> &x)
{
    // return x[0] * x[0] + x[1] * x[1] - 2 * x[0] - 4 * x[1] - 1; // F1

    // return 3 * x[0] * x[0] - 12 * x[0] + 2 * x[1] * x[1] + 16 * x[1] - 10; // F2

    // return x[0] * x[0] - 4 * x[0] * x[1] + 5 * x[1] * x[1] - 4 * x[1] + 3; // F3
    // ^ aceasta gaseste solutia doar cu backtracking = false (intr-un numar mare de pasi)

    return x[0] * x[0] * x[1] - 2 * x[0] * x[1] * x[1] + 3 * x[0] * x[1] + 4; // F4
}

// Gradient analitic al functiei de test
std::vector<double> gradientAnalitic(const std::vector<double> &x)
{
    // return {2 * x[0] - 2, 2 * x[1] - 4}; // F1
    // return {6 * x[0] - 12, 4 * x[1] + 16}; // F2
    // return {2 * x[0] - 4 * x[1], -4 * x[0] + 10 * x[1] - 4}; // F3
    return {2 * x[0] * x[1] - 2 * x[1] * x[1] + 3 * x[1], x[0] * x[0] - 4 * x[0] * x[1] + 3 * x[0]}; // F4
}

// Gradient numeric cu formula aproximativa data
std::vector<double> gradientNumeric(const std::vector<double> &x, std::function<double(const std::vector<double> &)> f)
{
    std::vector<double> gradient(n_dim);
    std::vector<double> x_mod = x; // Vectorul in care modificam doar valoarea indexului i

    for (int i = 0; i < n_dim; ++i)
    {
        x_mod[i] = x[i] + 2 * h;
        double F1 = f(x_mod);

        x_mod[i] = x[i] + h;
        double F2 = f(x_mod);

        x_mod[i] = x[i] - h;
        double F3 = f(x_mod);

        x_mod[i] = x[i] - 2 * h;
        double F4 = f(x_mod);

        x_mod[i] = x[i]; // Resetam pentru pasii urmatori

        gradient[i] = (-F1 + 8 * F2 - 8 * F3 + F4) / (12 * h);
    }

    return gradient;
}

// Norma unui vector
double norm(const std::vector<double> &v)
{
    double sum = 0;
    for (double val : v)
    {
        sum += val * val;
    }
    return std::sqrt(sum);
}

// Inmultirea unui vector cu un scalar
std::vector<double> inmultireScalar(const std::vector<double> &v, double s)
{
    std::vector<double> res(v.size());
    for (int i = 0; i < v.size(); ++i)
    {
        res[i] = v[i] * s;
    }
    return res;
}

// Scaderea a doi vectori
std::vector<double> scadereVectori(const std::vector<double> &a, const std::vector<double> &b)
{
    std::vector<double> res(a.size());
    for (int i = 0; i < a.size(); ++i)
    {
        res[i] = a[i] - b[i];
    }
    return res;
}

// Gradient descendent
std::vector<double> gradientDescendent(
    std::function<double(const std::vector<double> &)> f,
    std::function<std::vector<double>(const std::vector<double> &)> f_gradient)
{
    // Initializare random x0
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-5.0, 5.0); // Interval arbitrar (este strict punctul de start din spatiul solutiilor)
    // Scop: sa nu porneasca prea departe de punctele de minim

    std::vector<double> x(n_dim);
    for (int i = 1; i < n_dim; ++i)
    {
        x[i] = dis(gen); // Toate valorile x sunt luate aleatoriu la inceput, in jurul lui x0 (fiindca altfel nu pot aplica F(x))
    }

    int k = 0;
    while (k < kmax)
    {
        std::vector<double> gradient = f_gradient(x);

        // Calculul ratei de invatare:
        double eta = 0.001; // Metoda 1: Constanta
        if (backtracking)   // Metoda 2: Schimbam rata de invatare in functie de contextul local (backtracking line search) -> pentru functii convexe
        {
            eta = 1.0;
            double beta = 0.8;
            int p = 1;
            while (f(scadereVectori(x, inmultireScalar(gradient, eta))) > f(x) - (eta / 2) * norm(gradient) * norm(gradient) && p < 8)
            {
                eta *= beta;
                p++;
            }
        }

        std::vector<double> x_next = scadereVectori(x, inmultireScalar(gradient, eta)); // x_k+1 = x_k - eta_k * gradient_F(x_k)

        if (norm(scadereVectori(x_next, x)) <= epsilon) // Conform observatiei:
        // Consideram x_k0 aprox. egal cu x* cand diferenta dintre doua elemente succesive e suficient de mica
        {
            std::cout << "Punct gasit in " << k << " pasi \n";
            return x_next;
        }

        x = x_next; // x = x - eta * gradient
        k++;

        if (eta * norm(gradient) <= epsilon || eta * norm(gradient) >= 10e10)
        {
            std::cout << "Divergenta \n";
            return x;
        }
    }

    std::cout << "Depasit numarul maxim de pasi \n";
    return x;
}

// ---------------- MAIN ----------------
int main()
{
    auto sol1 = gradientDescendent(F, gradientAnalitic);
    std::cout << "Solutie aproximata (gradient analitic): (" << sol1[0] << ", " << sol1[1] << ")\n";

    auto sol2 = gradientDescendent(F, [&](const std::vector<double> &x)
                                   { return gradientNumeric(x, F); });
    std::cout << "Solutie aproximata (gradient numeric): (" << sol2[0] << ", " << sol2[1] << ")\n";

    return 0;
}