#include "Eigen/Dense"
#include <iostream>
#include <vector>

int main() {
    using namespace Eigen;
    using namespace std;

    // Example: dimensions p > n
    constexpr int p = 5;
    constexpr int n = 3;

    // Define a p x n matrix A
    Matrix<double, p, n> A;
    A << 1, 2, 3,
         4, 5, 6,
         7, 8, 9,
         2, 4, 6,
         3, 6, 9;

    // Define vector b of size p
    VectorXd b(p);
    b << 1, 2, 3, 4, 5;

    // Compute SVD
    JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);

    VectorXd singularValues = svd.singularValues();
    MatrixXd U = svd.matrixU();
    MatrixXd V = svd.matrixV();
    MatrixXd S = singularValues.asDiagonal();

    cout << "Singular values:\n" << singularValues << "\n\n";

    // Compute rank: number of non-zero singular values (with tolerance)
    double tol = numeric_limits<double>::epsilon() * max(A.rows(), A.cols()) * singularValues.array().abs()(0);
    int rank = (singularValues.array() > tol).count();

    cout << "Rank of A: " << rank << "\n\n";

    // Condition number: ratio of largest to smallest non-zero singular value
    double cond_number = singularValues(0) / singularValues(singularValues.size() - 1);
    cout << "Condition number of A: " << cond_number << "\n\n";

    // Compute pseudoinverse: AI = V * S_pinv * U.transpose()
    VectorXd singularValuesInv = singularValues;
    for (int i = 0; i < singularValues.size(); ++i) {
        if (singularValues(i) > tol)
            singularValuesInv(i) = 1.0 / singularValues(i);
        else
            singularValuesInv(i) = 0.0;
    }
    
    MatrixXd S_pinv = singularValuesInv.asDiagonal();
    MatrixXd A_pinv = V * S_pinv * U.transpose();

    cout << "Pseudoinverse of A (AI):\n" << A_pinv << "\n\n";

    // Compute xI = AI * b
    VectorXd xI = A_pinv * b;

    cout << "Solution xI = AI * b:\n" << xI << "\n\n";

    // Compute residual norm ||b - AxI||_2
    double residual_norm = (b - A * xI).norm();
    cout << "Residual norm ||b - AxI||_2: " << residual_norm << "\n";

    return 0;
}