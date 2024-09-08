#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

using namespace std;

// Function to create the tridiagonal matrix
vector<vector<double>> create_tridiagonal_matrix(int n) {
    vector<vector<double>> matrix(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; i++) {
        matrix[i][i] = 2.0004; // Main diagonal
        if (i > 0) {
            matrix[i][i-1] = -1; // Sub-diagonal
        }
        if (i < n-1) {
            matrix[i][i+1] = -1; // Super-diagonal
        }
    }

    matrix[n-1][n-2] = -2; // Adjust the last sub-diagonal element
    return matrix;
}

// Function to solve the system using the Gauss-Seidel method
vector<double> gauss_seidel(const vector<vector<double>>& A, const vector<double>& b, const vector<double>& x0, double tol = 0.0001, int max_iterations = 2498) {
    int n = b.size();
    vector<double> x = x0;
    vector<double> x_old(n);

    for (int iteration = 0; iteration < max_iterations; iteration++) {
        x_old = x;

        for (int i = 0; i < n; i++) {
            double sum1 = 0.0, sum2 = 0.0;

            // Sum for elements before i (lower triangular part)
            for (int j = 0; j < i; j++) {
                sum1 += A[i][j] * x[j];
            }

            // Sum for elements after i (upper triangular part)
            for (int j = i + 1; j < n; j++) {
                sum2 += A[i][j] * x_old[j];
            }

            // Update x[i]
            x[i] = (b[i] - sum1 - sum2) / A[i][i];
        }

        // Check for convergence
        double max_diff = 0.0;
        for (int i = 0; i < n; i++) {
            max_diff = max(max_diff, fabs(x[i] - x_old[i]));
        }

        if (max_diff < tol) {
            cout << "Converged in " << iteration + 1 << " iterations." << endl;
            return x;
        }
    }

    cout << "Maximum iterations reached without convergence." << endl;
    return x;
}

int main() {
    int n = 100;
    vector<vector<double>> A = create_tridiagonal_matrix(n);

    // Initial guess
    vector<double> x0(n, 1.0);

    // Right-hand side vector
    vector<double> b(n, 0.0);
    b[0] = 1.0;

    // Solve using Gauss-Seidel
    vector<double> solution = gauss_seidel(A, b, x0);

    // Output the solution
    cout << "Solution: ";
    for (double x : solution) {
        cout << x << " ";
    }
    cout << endl;

    return 0;
}
