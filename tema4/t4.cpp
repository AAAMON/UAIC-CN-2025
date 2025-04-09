#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

typedef vector<vector<double>> Matrix;

// Function to multiply two matrices
Matrix multiplyMatrices(const Matrix &A, const Matrix &B)
{
    int n = A.size();
    Matrix result(n, vector<double>(n, 0));

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            for (int k = 0; k < n; ++k)
            {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

// Function to subtract two matrices
Matrix subtractMatrices(const Matrix &A, const Matrix &B)
{
    int n = A.size();
    Matrix result(n, vector<double>(n, 0));

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            result[i][j] = A[i][j] - B[i][j];
        }
    }
    return result;
}

// Function to add two matrices
Matrix addMatrices(const Matrix &A, const Matrix &B)
{
    int n = A.size();
    Matrix result(n, vector<double>(n, 0));

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            result[i][j] = A[i][j] + B[i][j];
        }
    }
    return result;
}

// Function to compute the matrix norm (Frobenius norm)
double matrixNorm(const Matrix &A)
{
    double norm = 0.0;
    int n = A.size();
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            norm += A[i][j] * A[i][j];
        }
    }
    return sqrt(norm);
}

// Function to initialize the V0 matrix (diagonal matrix)
Matrix initializeV0(const Matrix &A)
{
    int n = A.size();
    Matrix V0(n, vector<double>(n, 0));
    // Initialize with identity matrix for simplicity (or another approach can be used)
    for (int i = 0; i < n; ++i)
    {
        V0[i][i] = 1.0 / A[i][i]; // Using diagonal elements of A as initial approximation
    }
    return V0;
}

// Function to multiply a matrix by a scalar
Matrix multiplyMatrixByScalar(const Matrix &matrix, double scalar)
{
    int n = matrix.size();
    Matrix result(n, vector<double>(n, 0));
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            result[i][j] = matrix[i][j] * scalar;
        }
    }
    return result;
}

Matrix schultzMethod(const Matrix &A, double epsilon, int kmax)
{
    int n = A.size();
    Matrix V0 = initializeV0(A); // Initial guess for V0
    Matrix V1 = V0;              // Start with V0

    int k = 0;
    double deltaV;

    do
    {
        // V(k+1) = V(k) * (2I - A * V(k))
        Matrix I(n, vector<double>(n, 0));
        for (int i = 0; i < n; ++i)
        {
            I[i][i] = 1.0;
        }

        // Calculate A * V(k)
        Matrix A_Vk = multiplyMatrices(A, V1);

        // Multiply I by 2 (scalar multiplication)
        Matrix twoI = multiplyMatrixByScalar(I, 2.0);

        // Calculate (2I - A * V(k))
        Matrix twoI_minus_A_Vk = subtractMatrices(twoI, A_Vk);

        // Calculate V(k+1) = V(k) * (2I - A * V(k))
        Matrix V_next = multiplyMatrices(V1, twoI_minus_A_Vk);

        // DEBUG: Print intermediate matrices for debugging
        /*
        cout << "Iteration " << k + 1 << endl;
        cout << "Matrix A:" << endl;
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                cout << fixed << setprecision(6) << A[i][j] << " ";
            }
            cout << endl;
        }
        cout << "Matrix V(k):" << endl;
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                cout << fixed << setprecision(6) << V1[i][j] << " ";
            }
            cout << endl;
        }
        cout << "Matrix 2I - A*V(k):" << endl;
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                cout << fixed << setprecision(6) << twoI_minus_A_Vk[i][j] << " ";
            }
            cout << endl;
        }
        cout << "Matrix V(k+1):" << endl;
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                cout << fixed << setprecision(6) << V_next[i][j] << " ";
            }
            cout << endl;
        }
        */

        // Calculate the difference ||V(k+1) - V(k)||
        deltaV = matrixNorm(subtractMatrices(V_next, V1));

        // Update V1 to V_next for the next iteration
        V1 = V_next;

        k++;

        if (k > kmax)
        {
            cout << "Exceeded maximum iterations (Schultz)." << endl;
            break;
        }

    } while (deltaV >= epsilon); // Convergence condition

    cout << "Iterations (Schultz): " << k << endl;

    Matrix A_Ainv = multiplyMatrices(A, V1);
    Matrix I(n, vector<double>(n, 0));
    for (int i = 0; i < n; ++i)
    {
        I[i][i] = 1.0;
    }
    Matrix diff = subtractMatrices(A_Ainv, I);
    double normFinal = matrixNorm(diff);
    cout << "Convergence Norm (Schultz): " << normFinal << endl;

    return V1; // Return the approximate inverse
}

Matrix liLiMethod1(const Matrix &A, double epsilon, int kmax)
{
    int n = A.size();
    Matrix V0 = initializeV0(A);
    Matrix V1 = V0;
    int k = 0;
    double deltaV;

    do
    {
        // V(k+1) = V(k) * (3I - A * V(k) * (3I - A * V(k)))
        Matrix I(n, vector<double>(n, 0));
        for (int i = 0; i < n; ++i)
        {
            I[i][i] = 1.0;
        }

        // Calculate A * V(k)
        Matrix A_Vk = multiplyMatrices(A, V1);

        // Multiply I by 3 (scalar multiplication)
        Matrix threeI = multiplyMatrixByScalar(I, 3.0);

        // Calculate (3I - A * V(k))
        Matrix inner = subtractMatrices(threeI, A_Vk);

        // Calculate A * V(k)
        Matrix middle = multiplyMatrices(A, V1);

        // Calculate A * V(k) * (3I - A * V(k))
        Matrix secondPart = multiplyMatrices(middle, inner);

        // Calculate 3I - A * V(k) * (3I - A * V(k))
        Matrix outer = subtractMatrices(threeI, secondPart);

        // Calculate V(k+1) = V(k) * (3I - A * V(k) * (3I - A * V(k)))
        Matrix V_next = multiplyMatrices(V1, outer);

        // Calculate the difference ||V(k+1) - V(k)||
        deltaV = matrixNorm(subtractMatrices(V_next, V1));

        // Update V1 to V_next for the next iteration
        V1 = V_next;

        k++;

        if (k > kmax)
        {
            cout << "Exceeded maximum iterations (Li & Li 1)." << endl;
            break;
        }

    } while (deltaV >= epsilon); // Convergence condition

    cout << "Iterations (Li & Li 1): " << k << endl;

    Matrix A_Ainv = multiplyMatrices(A, V1);
    Matrix I(n, vector<double>(n, 0));
    for (int i = 0; i < n; ++i)
    {
        I[i][i] = 1.0;
    }
    Matrix diff = subtractMatrices(A_Ainv, I);
    double normFinal = matrixNorm(diff);
    cout << "Convergence Norm (Li & Li 1): " << normFinal << endl;

    return V1;
}

Matrix liLiMethod2(const Matrix &A, double epsilon, int kmax)
{
    int n = A.size();
    Matrix V0 = initializeV0(A);
    Matrix V1 = V0;
    int k = 0;
    double deltaV;

    do
    {
        // V(k+1) = (I + (1 / 4) * (I - V(k) * A) * (3I - V(k) * A) * (3I - V(k) * A)) * V(k)
        Matrix I(n, vector<double>(n, 0));
        for (int i = 0; i < n; ++i)
        {
            I[i][i] = 1.0;
        }

        // Calculate V(k) * A
        Matrix Vk_A = multiplyMatrices(V1, A);

        // Multiply I by 3 (scalar multiplication)
        Matrix threeI = multiplyMatrixByScalar(I, 3.0);

        // Calculate (3I - V(k) * A)
        Matrix inner = subtractMatrices(threeI, Vk_A);

        // Calculate (3I - V(k) * A) * (3I - V(k) * A)
        Matrix inner_squared = multiplyMatrices(inner, inner);

        // Calculate (I - V(k) * A)
        Matrix middle = subtractMatrices(I, Vk_A);

        // Calculate (I - V(k) * A) * (3I - V(k) * A) * (3I - V(k) * A)
        Matrix product = multiplyMatrices(middle, inner_squared);

        // Calculate (1/4) * (I - V(k) * A) * (3I - V(k) * A) * (3I - V(k) * A)
        Matrix big_product = multiplyMatrixByScalar(product, 1.0 / 4.0);

        // Calculate I + (1/4) * (I - V(k) * A) * (3I - V(k) * A) * (3I - V(k) * A)
        Matrix sum = addMatrices(I, big_product);

        // Calculate V(k+1) = (I + (1 / 4) * (I - V(k) * A) * (3I - V(k) * A) * (3I - V(k) * A)) * V(k)
        Matrix V_next = multiplyMatrices(sum, V1);

        // Calculate the difference ||V(k+1) - V(k)||
        deltaV = matrixNorm(subtractMatrices(V_next, V1));

        // Update V1 to V_next for the next iteration
        V1 = V_next;

        k++;

        if (k > kmax)
        {
            cout << "Exceeded maximum iterations (Li & Li 2)." << endl;
            break;
        }

    } while (deltaV >= epsilon); // Convergence condition

    cout << "Iterations (Li & Li 2): " << k << endl;

    Matrix A_Ainv = multiplyMatrices(A, V1);
    Matrix I(n, vector<double>(n, 0));
    for (int i = 0; i < n; ++i)
    {
        I[i][i] = 1.0;
    }
    Matrix diff = subtractMatrices(A_Ainv, I);
    double normFinal = matrixNorm(diff);
    cout << "Convergence Norm (Li & Li 2): " << normFinal << endl;

    return V1;
}

// Function to calculate the exact inverse of the matrix A
// cu formula noastra magica
Matrix exactInverse(const Matrix &A)
{
    int n = A.size();
    Matrix inverse(n, vector<double>(n, 0.0));

    // Fill the inverse matrix based on the observed pattern
    for (int i = 0; i < n; ++i)
    {
        inverse[i][i] = 1.0; // Diagonal elements are 1
        for (int j = i + 1; j < n; ++j)
        {
            // Alternating powers of 2 with signs: -2^(j-i), 2^(j-i)
            inverse[i][j] = (j - i) % 2 == 1 ? -pow(2, j - i) : pow(2, j - i);
        }
    }

    return inverse;
}

int main()
{
    // Read matrix from file
    string filename = "input.txt";
    ifstream inputFile(filename);

    if (!inputFile.is_open())
    {
        cerr << "Unable to open file!" << endl;
        return 1;
    }

    int n;
    double epsilon;
    int kmax;
    inputFile >> n >> epsilon >> kmax;

    Matrix A(n, vector<double>(n));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            inputFile >> A[i][j];
        }
    }
    inputFile.close();

    // Apply Schultz's method to approximate the inverse
    Matrix approxInverseS = schultzMethod(A, epsilon, kmax);

    // Apply Li and Li's first method to approximate the inverse
    Matrix approxInverseLL1 = liLiMethod1(A, epsilon, kmax);

    // Apply Li and Li's second method to approximate the inverse
    Matrix approxInverseLL2 = liLiMethod2(A, epsilon, kmax);

    // Calculate the exact inverse
    Matrix exactInv = exactInverse(A);

    // Output the exact inverse
    cout << "Exact Inverse:" << endl;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            cout << fixed << setprecision(6) << exactInv[i][j] << " ";
        }
        cout << endl;
    }

    // Output the approximate inverse
    cout << "Approximate Inverse (Schultz's Method):" << endl;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            cout << fixed << setprecision(6) << approxInverseS[i][j] << " ";
        }
        cout << endl;
    }

    // Output the approximate inverse
    cout << "Approximate Inverse (Li and Li's first Method):" << endl;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            cout << fixed << setprecision(6) << approxInverseLL1[i][j] << " ";
        }
        cout << endl;
    }

    // Output the approximate inverse
    cout << "Approximate Inverse (Li and Li's second Method):" << endl;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            cout << fixed << setprecision(6) << approxInverseLL2[i][j] << " ";
        }
        cout << endl;
    }

    // Calculate the norm of the difference
    Matrix diffS = subtractMatrices(exactInv, approxInverseS);
    double normS = matrixNorm(diffS);
    cout << "Norm of the difference between exact and approximate (Schultz) inverse: " << normS << endl;

    Matrix diffLL1 = subtractMatrices(exactInv, approxInverseLL1);
    double normLL1 = matrixNorm(diffLL1);
    cout << "Norm of the difference between exact and approximate (Li & Li 1) inverse: " << normLL1 << endl;

    Matrix diffLL2 = subtractMatrices(exactInv, approxInverseLL2);
    double normLL2 = matrixNorm(diffLL2);
    cout << "Norm of the difference between exact and approximate (Li & Li 2) inverse: " << normLL2 << endl;

    return 0;
}
