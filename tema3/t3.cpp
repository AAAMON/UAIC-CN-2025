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

// Function to read vector b from file
vector<double> readVector(const string& fileName) {
    ifstream file(fileName);
    if (!file.is_open()) {
        throw runtime_error("Could not open file: " + fileName);
    }

    vector<double> b;
    int n;  // Dimension of the vector
    file >> n;
    b.resize(n);
    
    for (int i = 0; i < n; ++i) {
        file >> b[i];
    }

    return b;
}

//////////////////////////////////////////////////////////////////////////
// 2 /////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

// Compressed Sparse Row (CSR) representation for sparse matrices
struct SparseMatrixCSR {
    int n;  // Dimension of the matrix
    vector<double> values;         // Non-zero values
    vector<int> colIndices;       // Column indices of the non-zero elements
    vector<int> rowPointers;      // Row pointers, indicating start of each row in values

    // Function to check if diagonal elements are non-zero (similar to the previous method)
    bool checkDiagonal() const {
        for (int i = 0; i < n; ++i) {
            bool found = false;
            for (int j = rowPointers[i]; j < rowPointers[i + 1]; ++j) {
                if (colIndices[j] == i) {
                    if (values[j] == 0.0) {
                        cout << "Diagonal element at position " << i << " is zero\n";
                        return true;
                    }
                    found = true;
                    break;
                }
            }
            if (!found) {
                cout << "Diagonal element at position " << i << " is missing\n";
                return false;
            }
        }
        return false;
    }
};

// Function to read a matrix in CSR format from file
SparseMatrixCSR readMatrixCSR(const string& fileName) {
    ifstream file(fileName);
    if (!file.is_open()) {
        throw runtime_error("Could not open file: " + fileName);
    }

    SparseMatrixCSR matrix;
    int n;  // Dimension of the matrix
    file >> n;
    matrix.n = n;

    // Initialize the rowPointers array
    matrix.rowPointers.resize(n + 1, 0);

    // Temporary storage for non-zero values and their positions
    vector<double> values;
    vector<int> colIndices;

    string line;
    getline(file, line);  // Skip the newline character after n

    // Read the non-zero elements (value, row, col) in the file
    while (getline(file, line)) {
        stringstream ss(line);
        double value;
        int row, col;
        char comma1, comma2;
        if (ss >> value >> comma1 >> row >> comma2 >> col) {

            // Store the non-zero elements in the CSR arrays
            values.push_back(value);
            colIndices.push_back(col);
            matrix.rowPointers[row + 1]++;
        }
    }

    // Accumulate the rowPointers (it stores the starting index of each row)
    for (int i = 1; i <= n; ++i) {
        matrix.rowPointers[i] += matrix.rowPointers[i - 1];
    }

    // Copy the non-zero values and column indices into the matrix
    matrix.values = move(values);
    matrix.colIndices = move(colIndices);

    return matrix;
}

//////////////////////////////////////////////////////////////////////////
// 3 /////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

// Gauss-Seidel Method for solving Ax = b
void gaussSeidel(const SparseMatrix& A, vector<double>& x, const vector<double>& b, int maxIter = 1000, double tol = 1e-6) {
    int n = A.n;
    vector<double> x_old(n, 0.0); // Initial guess (starting with 0 for all variables)
    
    for (int iter = 0; iter < maxIter; ++iter) {
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            
            // Calculate the sum of A_ij * x_j for all j, excluding i
            for (const auto& element : A.rows[i].elements) {
                int j = element.first;
                double value = element.second;
                if (j != i) {
                    sum += value * x[j];
                }
            }
            
            // Calculate the new value for x[i] using the Gauss-Seidel formula
            x[i] = (b[i] - sum) / A.diagonal[i];
        }
        
        // Check for convergence by comparing the difference between new and old x
        double maxDiff = 0.0;
        for (int i = 0; i < n; ++i) {
            maxDiff = max(maxDiff, abs(x[i] - x_old[i]));
        }

        if (maxDiff < tol) {
            cout << "Converged after " << iter + 1 << " iterations.\n";
            return;
        }

        // Update x_old with the current x values for the next iteration
        x_old = x;
    }

    cout << "Maximum iterations reached. Solution may not have converged.\n";
}


int main() {
    try {
        for (int i = 1; i < 6; ++i)
        {
            // reprezentarea 1
            SparseMatrix A = readMatrix("data/a_" + to_string(i) + ".txt");
            vector<double> B = readVector("data/b_" + to_string(i) + ".txt");

            if (A.checkDiagonal()) {
                cout << "Matrix A" << i << " has zero on the diagonal!" << endl;
            } else {
                cout << "Matrix A" << i << "  has no zero on the diagonal." << endl;
            }

            // reprezentarea 2
            SparseMatrixCSR _A = readMatrixCSR("data/a_" + to_string(i) + ".txt");

            if (_A.checkDiagonal()) {
                cout << "Matrix _A" << i << " has zero on the diagonal!" << endl;
            } else {
                cout << "Matrix _A" << i << "  has no zero on the diagonal." << endl;
            }

            // calculam caca cu Gauss-Seidel
            vector<double> x(A.n, 0.0); // Initial solution vector (starts with 0)
            gaussSeidel(A, x, B); // Solve using Gauss-Seidel
            // Output the result
            cout << "Solution x = [";
            for (double val : x) {
                cout << val << " ";
            }
            cout << "]\n";

        }
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
    }

    return 0;
}
