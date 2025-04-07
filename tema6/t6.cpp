#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include <random>
#include <cmath>

using namespace std;
using namespace Eigen;

// Example function f(x)
double f(double x)
{
    return pow(x, 4) - 12 * pow(x, 3) + 30 * pow(x, 2) + 12;
}

// Generate n+1 sorted distinct points in (x0, xn)
vector<double> generate_x_values(double x0, double xn, int n)
{
    vector<double> x_values(n + 1);
    x_values[0] = x0;
    x_values[n] = xn;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(x0, xn);

    for (int i = 1; i < n; ++i)
    {
        double val;
        do
        {
            val = dist(gen);
        } while (val <= x_values[i - 1]);
        x_values[i] = val;
    }
    return x_values;
}

// Construct the Vandermonde matrix for least squares
MatrixXd construct_vandermonde(const vector<double> &x_values, int m)
{
    int n = x_values.size();
    MatrixXd V(n, m + 1);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j <= m; ++j)
        {
            V(i, j) = pow(x_values[i], j);
        }
    }
    return V;
}

// Compute polynomial coefficients using least squares
VectorXd least_squares_fit(const vector<double> &x_values, const vector<double> &y_values, int m)
{
    MatrixXd V = construct_vandermonde(x_values, m);
    VectorXd y = Map<const VectorXd>(y_values.data(), y_values.size());
    VectorXd coeffs = (V.transpose() * V).ldlt().solve(V.transpose() * y);
    return coeffs;
}

// Evaluate polynomial at x using Horner's method
double horner_eval(const VectorXd &coeffs, double x)
{
    double result = 0.0;
    for (int i = coeffs.size() - 1; i >= 0; --i)
    {
        result = result * x + coeffs[i];
    }
    return result;
}

// Interpolare trigonometrica
void interpolateTrig(int m, double x_bar)
{
    int n = 2 * m + 1; // numar impar de puncte
    vector<double> x(n);
    vector<double> y(n);

    // Generam puncte echidistante in [0, 2pi)
    for (int i = 0; i < n; ++i)
    {
        x[i] = 2 * M_PI * i / n;
        y[i] = f(x[i]);
    }

    // Matricea T va avea dimensiunea n x (2m+1)
    MatrixXd T(n, 2 * m + 1);

    for (int i = 0; i < n; ++i)
    {
        T(i, 0) = 1; // a0 (cos(0x) = 1)
        for (int k = 1; k <= m; ++k)
        {
            T(i, 2 * k - 1) = cos(k * x[i]); // coloana a_k
            T(i, 2 * k) = sin(k * x[i]);     // coloana b_k
        }
    }

    // Vectorul coloana cu valorile functiei
    VectorXd Y(n);
    for (int i = 0; i < n; ++i)
    {
        Y(i) = y[i];
    }

    // Rezolva sistemul T * X = Y cu metoda QR
    VectorXd coeffs = T.colPivHouseholderQr().solve(Y);

    // Evaluam T_n(x_bar)
    double Tn_x_bar = coeffs(0); // a0
    for (int k = 1; k <= m; ++k)
    {
        Tn_x_bar += coeffs(2 * k - 1) * cos(k * x_bar); // a_k * cos(kx)
        Tn_x_bar += coeffs(2 * k) * sin(k * x_bar);     // b_k * sin(kx)
    }

    // Calculam eroarea punctuala si eroarea totala
    double err_point = abs(f(x_bar) - Tn_x_bar);
    double total_err = 0.0;
    for (int i = 0; i < n; ++i)
    {
        double Ti = coeffs(0);
        for (int k = 1; k <= m; ++k)
        {
            Ti += coeffs(2 * k - 1) * cos(k * x[i]);
            Ti += coeffs(2 * k) * sin(k * x[i]);
        }
        total_err += abs(y[i] - Ti);
    }

    cout << "\n[Interpolare Trigonometrica]" << endl;
    cout << "T_n(x_bar) = " << Tn_x_bar << endl;
    cout << "|T_n(x_bar) - f(x_bar)| = " << err_point << endl;
    cout << "Total error: " << total_err << endl;
}


// Enter x0, xn (x0 < xn): 0
// 6.28
// Enter number of points (n): 20
// Enter polynomial degree (m < 6): 4
// Enter x_bar (point to approximate): 1.5





int main()
{
    int n, m;
    double x0, xn, x_bar;

    cout << "Enter x0, xn (x0 < xn): ";
    cin >> x0 >> xn;
    cout << "Enter number of points (n): ";
    cin >> n;
    cout << "Enter polynomial degree (m < 6): ";
    cin >> m;
    cout << "Enter x_bar (point to approximate): ";
    cin >> x_bar;

    vector<double> x_values = generate_x_values(x0, xn, n);
    vector<double> y_values(n + 1);
    for (int i = 0; i <= n; ++i)
    {
        y_values[i] = f(x_values[i]);
    }

    VectorXd coeffs = least_squares_fit(x_values, y_values, m);
    double Pm_x_bar = horner_eval(coeffs, x_bar);

    double error_approx = abs(Pm_x_bar - f(x_bar));
    double total_error = 0.0;
    for (int i = 0; i <= n; ++i)
    {
        total_error += abs(horner_eval(coeffs, x_values[i]) - y_values[i]);
    }

    cout << "P_m(x_bar) = " << Pm_x_bar << endl;
    cout << "|P_m(x_bar) - f(x_bar)| = " << error_approx << endl;
    cout << "Total error: " << total_error << endl;

    interpolateTrig(m, x_bar);
    return 0;
}