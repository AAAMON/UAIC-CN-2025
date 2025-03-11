#include <iostream>
#include <vector>
#include <cmath>
//#include <Eigen/Dense>

using namespace std;

bool abortt = 0;


bool isZero(double x, double epsilon) {
    if (abs(x) <= epsilon)
    {
        cerr << "Error! " << x << " is almost 0 (zero), cannot divide by " << x << '\n';
        return 1;
    }
    return 0;
}

void descompunereLU(vector<vector<double>>& A, vector<double>& dU, int n, double epsilon) {
    
    // PAS 1 - Descompunere LU

    // Algoritmul Crout
    // la pasul p calculam:
    // - elementele de pe COLOANA p pentru L
    // - elementele de pe LINIA p pentru U

    for (int p = 0; p < n; ++p) {

        // L
        // Pentru L trebuie sa calculam si ce e pe diagonala, motiv pentru care
        // i incepe cu valoarea p in loc de p-1

        for (int i = p; i < n; ++i) {
            // scadem ccturi
            // don't worry about it ok it works
            for (int z = 0; z < p; ++z) {
                A[i][p] -= A[i][z]*A[z][p];
            }
            // TODO: la toate impartirile sa verificam ca nu e 0 folosing epsilon
            if (!isZero(dU[p], epsilon))
                A[i][p] /= dU[p];
            else {
                abortt = true;
                return;
            }
        }

        // U
        // Pentru U deja avem ce are el pe diagonala in dU
        for (int i = p+1; i < n; ++i) {
            // scadem ccturi
            for (int z = 0; z < p; ++z) {
                A[p][i] -= A[p][z]*A[z][i];
            }
            // TODO: la toate impartirile sa verificam ca nu e 0 folosing epsilon
            if (!isZero(A[p][p], epsilon))
                A[p][i] /= A[p][p];
            else {
                abortt = true;
                return;
            }
        }
    }

    // Printam descompunerea
    cout << "Matrix:" << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
}

double calculDeterminant(vector<vector<double>>& A, vector<double>& dU, int n, double epsilon) {
    
    // Pas 2 - Calcul det A

    // det A = det L * det U

    double detA = 1;
    for (int i = 0; i < n; ++i) {
        detA *= A[i][i] * dU[i];
    }
    
    cout << "Det A = " << detA << '\n';

    return detA;
}

void run(vector<vector<double>>& A, vector<double>& dU, int n, double epsilon)
{
    descompunereLU(A, dU, n, epsilon);
    if (abortt == true) {
        return;
    }
    double detA = calculDeterminant(A, dU, n, epsilon);
}


int main () 
{
    // Dimensiunea sistemului
    int n = 3;

    // Precizia
    double epsilon = 1e-10;

    // Matricea A
    vector<vector<double>> A = {
        {2.5, 2, 2},
        {-5, -2, -3},
        {5, 6, 6.5}
    };

    // Vectorul dU
    vector<double> dU = {1, 1, 1}; 

    run(A, dU, n, epsilon);



    return 0;
}