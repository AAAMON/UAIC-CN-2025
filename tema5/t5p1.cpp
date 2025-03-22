#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <stdexcept>
#include <sstream>
#include <string>
#include <map>
#include <cmath>  // For abs() to check convergence
#include <limits>  // For the max value used in checking convergence
#include <algorithm>  // For find_if
#include <random>     // For random number generation

using namespace std;

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
                return true;  // Return the index of the zero element
            }
        }
        return false;  // Return -1 if no zero element is found
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
                    return false; // Matrix is not symmetric
                }
            }
        }
        return true; // Matrix is symmetric
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

// Function to generate a random sparse symmetric matrix
SparseMatrix generateRandomSparseSymmetricMatrix(int size, double sparsity = 0.1, double valueRangeMin = 0.1, double valueRangeMax = 10.0) {
    SparseMatrix matrix;
    matrix.n = size;
    matrix.diagonal.resize(size, 0.0);
    matrix.rows.resize(size);

    // Random number generation setup
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(valueRangeMin, valueRangeMax);
    bernoulli_distribution sparsityDis(sparsity);

    // Fill the matrix
    for (int i = 0; i < size; ++i) {
        // Diagonal element
        matrix.diagonal[i] = dis(gen);

        for (int j = i + 1; j < size; ++j) {
            if (sparsityDis(gen)) {
                double value = dis(gen);
                matrix.rows[i].elements.push_back({j, value});
                matrix.rows[j].elements.push_back({i, value}); // Ensure symmetry
            }
        }
    }
    return matrix;
}

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

// Main function for testing
int main() {
    try {
        // Generate random sparse symmetric matrix and save to file
        int n = 500;  // Size of the matrix
        SparseMatrix randomMatrix = generateRandomSparseSymmetricMatrix(n);

        // Write the matrix to a file (for demonstration purposes)
        ofstream outFile("outputp1.txt");
        outFile << n << endl;
        for (int i = 0; i < n; ++i) {
            outFile << randomMatrix.diagonal[i] << " , " << i << " , " << i << endl;  // Diagonal elements
            for (const auto& elem : randomMatrix.rows[i].elements) {
                outFile << elem.second << " , " << i << " , " << elem.first << endl;
            }
        }

        // Read matrix from file and check symmetry
        SparseMatrix A = readMatrix("outputp1.txt");

        if (A.isSymmetric()) {
            cout << "Matrix is symmetric." << endl;
        } else {
            cout << "Matrix is not symmetric." << endl;
        }

    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
    }

    return 0;
}
