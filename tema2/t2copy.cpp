#include <iostream>
#include <fstream>
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
    int n = 100;
    double epsilon = 1e-10;

    // Randomly seed the random number generator
    srand(time(0));

    // Create a 100x100 sparse matrix (mostly zeros, with random values on the diagonal and some random off-diagonal values)
    vector<vector<double>> A(n, vector<double>(n, 0));  // Initialize a 100x100 matrix with zeros

    // Fill the diagonal with non-zero values
    for (int i = 0; i < n; i++) {
        A[i][i] = (rand() % 100 + 1);  // Assign a random non-zero value to the diagonal
    }

    // Optionally fill some off-diagonal values to make the matrix sparse
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j && rand() % 5 == 0) {  // Add non-zero elements with a certain probability
                A[i][j] = rand() % 10 + 1;  // Random off-diagonal value
            }
        }
    }

    // Create vector b with random values
    vector<double> b(n);
    for (int i = 0; i < n; i++) {
        b[i] = rand() % 100 + 1;  // Random values for vector b
    }

    // Create vector dU with ones (or you can customize it)
    //vector<double> dU(n, 1);  // Initialize vector dU with 1's (you can adjust this if necessary)
    vector<double> dU(n);
    for (int i = 0; i < n; i++) {
        dU[i] = rand() % 100 + 1;  // Random values for vector b
    }

    // Save to file "input.txt"
    ofstream outputFile("output.txt");
    if (!outputFile) {
        cerr << "Failed to open the file!" << endl;
        return 1;
    }

    // Write n and epsilon
    outputFile << n << endl;
    outputFile << epsilon << endl;

    // Write matrix A to file
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            outputFile << A[i][j] << " ";
        }
        outputFile << endl;
    }

    // Write vector b to file
    for (int i = 0; i < n; i++) {
        outputFile << b[i] << " ";
    }
    outputFile << endl;

    // Write vector dU to file
    for (int i = 0; i < n; i++) {
        outputFile << dU[i] << " ";
    }
    outputFile << endl;

    // Close the file
    outputFile.close();


    run(A, dU, b, n, epsilon);

    return 0;
}