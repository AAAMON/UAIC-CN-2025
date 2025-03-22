#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <stdexcept>
#include <sstream>
#include <string>
#include <map>
#include <cmath> // For abs() to check convergence
#include <limits> // For the max value used in checking convergence
#include <algorithm> // For find_if


using namespace std;

//////////////////////////////////////////////////////////////////////////
// 1 /////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

// Structure for storing sparse row
struct SparseRow {
    vector<pair<int, double>> elements;  // (column index, value) for non-zero elements
};

// Sparse matrix structure with diagonal vector and sparse row vectors
struct SparseMatrix {
    int n;  // Dimension of the matrix
    vector<double> diagonal;  // Diagonal elements of the matrix
    vector<SparseRow> rows;  // Sparse rows (excluding diagonal elements)

    // Function to check if diagonal elements are non-zero
    bool checkDiagonal() const {
        for (int i = 0; i < diagonal.size(); ++i) {
            if (diagonal[i] == 0.0) {
                cout << "On i: " << i << '\n';
                return 1;  // Return the index of the zero element
            }
        }
        return 0;  // Return -1 if no zero element is found
    }

    // Function to check if the matrix is symmetric
    bool isSymmetric() const {
        for (int i = 0; i < n; ++i) {
            // Check if off-diagonal elements (i, j) are equal to (j, i)
            for (const auto& elem : rows[i].elements) {
                int col = elem.first;
                double value = elem.second;
                auto it = find_if(rows[col].elements.begin(), rows[col].elements.end(),
                                  [i](const pair<int, double>& p) { return p.first == i; });
                if (it == rows[col].elements.end() || abs(it->second - value) > 1e-6) {
                    return false; // Matricea nu este simetrică
                }
            }
        }
        return true; // Matricea este simetrică
    }

    // Function to multiply matrix with a vector
    vector<double> multiplyWithVector(const vector<double>& vec) const {
        vector<double> result(n, 0.0);
        for (int i = 0; i < n; ++i) {
            result[i] = diagonal[i] * vec[i];  // Multiply diagonal
            for (const auto& elem : rows[i].elements) {
                result[i] += elem.second * vec[elem.first];  // Multiply non-diagonal elements
            }
        }
        return result;
    }
};

// Function to read matrix A from file (updated for your format: value, row, column)
SparseMatrix readMatrix(const string& fileName) {
    ifstream file(fileName);
    if (!file.is_open()) {
        throw runtime_error("Could not open file: " + fileName);
    }

    SparseMatrix matrix;
    int n;  // Dimension of the matrix
    file >> n;
    matrix.n = n;

    matrix.diagonal.resize(n, 0.0);
    matrix.rows.resize(n);

    string line;
    getline(file, line);  // Read the rest of the first line (to handle newline after n)

    while (getline(file, line)) {  // Read file line by line
        stringstream ss(line);  // Use stringstream to process the line
        double value;
        int row, col;

        // Read the value, row, and col, expecting the format "value , row , col"
        char comma1, comma2;
        if (ss >> value >> comma1 >> row >> comma2 >> col) {
            // Handle malformed lines
            if (comma1 != ',' || comma2 != ',' || row < 0 || col < 0) {
                continue;
            }
            
            if (row == col) {
                matrix.diagonal[row] = value;  // Diagonal elements are stored in the diagonal vector
            } else {
                matrix.rows[row].elements.push_back({col, value});  // Non-diagonal elements
            }
        } else {
            // Handle malformed line (skip it)
            continue;
        }
    }

    return matrix;
}


// Function to apply the power method to approximate the largest eigenvalue and eigenvector
pair<double, vector<double>> powerMethod(const SparseMatrix& matrix, int maxIterations = 1000, double tol = 1e-6) {
    vector<double> u(matrix.n, 1.0);  // Initial guess for eigenvector
    double lambda = 0.0;
    double prevLambda = 0.0;

    for (int iter = 0; iter < maxIterations; ++iter) {
        vector<double> Au = matrix.multiplyWithVector(u);
        
        // Find the largest value in Au to normalize the vector
        double norm = 0.0;
        for (double val : Au) {
            norm += val * val;
        }
        norm = sqrt(norm);
        
        // Normalize the vector
        for (int i = 0; i < matrix.n; ++i) {
            u[i] = Au[i] / norm;
        }
        
        // Approximate the largest eigenvalue
        lambda = 0.0;
        for (int i = 0; i < matrix.n; ++i) {
            lambda += u[i] * Au[i];
        }

        // Check for convergence
        if (abs(lambda - prevLambda) < tol) {
            break;
        }

        prevLambda = lambda;
    }

    return {lambda, u};
}

// Function to calculate the norm of a vector
double calculateNorm(const vector<double>& vec1, const vector<double>& vec2) {
    double norm = 0.0;
    for (int i = 0; i < vec1.size(); ++i) {
        norm += (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
    }
    return sqrt(norm);
}

int main() {
    try {
        for (int i = 1; i < 8; ++i) {
            // Read the matrix from file
            SparseMatrix A = readMatrix("data/a_" + to_string(i) + ".txt");

            // Check if the matrix is symmetric
            if (A.isSymmetric()) {
                cout << "Matrix A" << i << " is symmetric." << endl;
            } else {
                cout << "Matrix A" << i << " is not symmetric." << endl;
                continue;
            }

            // Apply power method to approximate the largest eigenvalue and eigenvector
            auto [lambdaMax, uMax] = powerMethod(A);

            // Calculate the norm |A * u_max - lambda_max * u_max|
            vector<double> Au = A.multiplyWithVector(uMax);
            for (int i = 0; i < uMax.size(); ++i) {
                Au[i] -= lambdaMax * uMax[i];
            }
            double norm = calculateNorm(Au, vector<double>(uMax.size(), 0.0));

            cout << "Matrix A" << i << " - Largest eigenvalue: " << lambdaMax << endl;
            cout << "Norm of |A * u_max - lambda_max * u_max|: " << norm << endl;
        }
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
    }

    return 0;
}
