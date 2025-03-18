#include <iostream>
#include <vector>
#include <cmath>
#include "Eigen/Dense"

using namespace std;

bool abortt = 0;

bool isZero(double x, double epsilon)
{
    if (abs(x) <= epsilon)
    {
        cerr << "Error! " << x << " is almost 0 (zero), cannot divide by " << x << '\n';
        return 1;
    }
    return 0;
}

void descompunereLU(vector<vector<double>> &A, vector<double> &dU, int n, double epsilon)
{

    // PAS 1 - Descompunere LU

    // Algoritmul Crout
    // la pasul p calculam:
    // - elementele de pe COLOANA p pentru L
    // - elementele de pe LINIA p pentru U

    for (int p = 0; p < n; ++p)
    {

        // L
        // Pentru L trebuie sa calculam si ce e pe diagonala, motiv pentru care
        // i incepe cu valoarea p in loc de p-1

        for (int i = p; i < n; ++i)
        {
            // scadem ccturi
            // don't worry about it ok it works
            for (int z = 0; z < p; ++z)
            {
                A[i][p] -= A[i][z] * A[z][p];
            }
            // TODO: la toate impartirile sa verificam ca nu e 0 folosing epsilon
            if (!isZero(dU[p], epsilon))
                A[i][p] /= dU[p];
            else
            {
                abortt = true;
                return;
            }
        }

        // U
        // Pentru U deja avem ce are el pe diagonala in dU
        for (int i = p + 1; i < n; ++i)
        {
            // scadem ccturi
            for (int z = 0; z < p; ++z)
            {
                A[p][i] -= A[p][z] * A[z][i];
            }
            // TODO: la toate impartirile sa verificam ca nu e 0 folosing epsilon
            if (!isZero(A[p][p], epsilon))
                A[p][i] /= A[p][p];
            else
            {
                abortt = true;
                return;
            }
        }
    }

    // Printam descompunerea
    cout << "Matrix:" << endl;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
}

double calculDeterminant(vector<vector<double>> &A, vector<double> &dU, int n, double epsilon)
{

    // Pas 2 - Calcul det A

    // det A = det L * det U

    double detA = 1;
    for (int i = 0; i < n; ++i)
    {
        detA *= A[i][i] * dU[i];
    }

    cout << "Det A = " << detA << '\n';

    return detA;
}

vector<double> substitutieDirecta(vector<vector<double>> &A, vector<double> b, int n, double epsilon) // pentru matrici inferior triunghiulare (Ly=b)
{
    vector<double> y(n, 0.0);

    for (int i = 0; i < n; i++)
    {
        double sum = 0.0;
        for (int j = 0; j < i; j++)
        {
            sum += A[i][j] * y[j];
        }
        if (!isZero(A[i][i], epsilon))
        {
            y[i] = (b[i] - sum) / A[i][i];
        }
        else
        {
            abortt = true;
            return y; // can't have empty return statement, but it doesnt matter anyway
        }
    }

    return y;
}

vector<double> substitutieInversa(vector<vector<double>> &A, vector<double> dU, vector<double> b, int n, double epsilon) // pentru matrici superior triunghiulare (Ux=y)
{
    vector<double> y = substitutieDirecta(A, b, n, epsilon);
    vector<double> x(n, 0.0);

    for (int i = n - 1; i >= 0; i--)
    {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++)
        {
            sum += A[i][j] * x[j];
        }
        if (!isZero(dU[i], epsilon)) // diagonala lui U se afla in dU. diagonala lui L a ramas in A
        {
            x[i] = (y[i] - sum) / dU[i];
        }
        else
        {
            abortt = true;
            return x;
        }
    }

    cout << "xLU = (";
    for (int i = 0; i < n - 1; i++)
    {
        cout << x[i] << ", ";
    }
    cout << x[n - 1] << ")\n";

    return x;
}

void calculNorma(const vector<vector<double>> &Ainit, vector<double> &xLU, vector<double> &b, int n)
{
    double sum = 0.0;
    for (int i = 0; i < n; i++)
    {
        double yi = 0.0;
        for (int j = 0; j < n; j++)
        {
            yi += Ainit[i][j] * xLU[j];
        }
        double zi = yi - b[i];
        // Dacă yi = b[i] pentru fiecare i, atunci zi = 0 și suma rămâne 0, ceea ce duce la norma 0
        // (xLU este o aproximare a solutiei, dar pentru dimensiuni mici rezultatul este exact)
        sum += zi * zi;
    }

    double norm = sqrt(sum);

    cout << "||Ainit * xLU - b||2 = " << norm << "\n";
}

void calculEigen(const vector<vector<double>> &Ainit, vector<double> &b, vector<double> xLU, int n)
{
    Eigen::MatrixXd Alib(n, n); // matrice de dimensiune arbitrara (n, n), cu valori double
    Eigen::VectorXd blib(n);

    for (int i = 0; i < n; i++)
    {
        blib(i) = b[i];
        for (int j = 0; j < n; j++)
        {
            Alib(i, j) = Ainit[i][j];
        }
    }

    Eigen::VectorXd xlib = Alib.colPivHouseholderQr().solve(blib); // solutia pentru A*x = b
    cout << "xlib = " << xlib.transpose() << "\n";
    Eigen::MatrixXd Alib_inv = Alib.inverse();
    cout << "A^(-1):\n"
         << Alib_inv << "\n";

    Eigen::VectorXd xLUlib(n);
    for (int i = 0; i < n; i++)
    {
        xLUlib(i) = xLU[i];
    }

    double norm1 = (xLUlib - xlib).norm();
    cout << "||xLU - xlib||2 = " << norm1 << "\n";

    Eigen::VectorXd xalt = Alib_inv * blib;
    double norm2 = (xLUlib - xalt).norm();
    cout << "||xLU - Alib_inv * b||2 = " << norm2 << "\n";
}

void run(vector<vector<double>> &A, vector<double> &dU, vector<double> &b, int n, double epsilon)
{
    const vector<vector<double>> Ainit = A;

    descompunereLU(A, dU, n, epsilon);
    if (abortt == true)
    {
        return;
    }

    double detA = calculDeterminant(A, dU, n, epsilon);

    vector<double> xLU = substitutieInversa(A, dU, b, n, epsilon);
    if (abortt == true)
    {
        return;
    }

    calculNorma(Ainit, xLU, b, n);

    calculEigen(Ainit, b, xLU, n);
}

int main()
{
    // Dimensiunea sistemului
    int n = 3;
    // std::cin >> n;

    // Precizia
    double epsilon = 1e-10;
    // std::cin >> epsilon;

    // Matricea A
    vector<vector<double>> A = {
        {2.5, 2, 2},
        {-5, -2, -3},
        {5, 6, 6.5}};

    // Vectorul termenilor liberi b
    vector<double> b = {2, -6, 2};

    // Vectorul dU
    vector<double> dU = {1, 1, 1};

    run(A, dU, b, n, epsilon);

    return 0;
}