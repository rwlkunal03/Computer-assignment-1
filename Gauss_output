#include <iostream>
#include <vector>

using namespace std;

// Thomas algorithm function for solving tridiagonal systems
vector<double> TDMAsolver(const vector<double>& a, const vector<double>& b, const vector<double>& c, const vector<double>& d) {
    int n = d.size();
    vector<double> ac = a, bc = b, cc = c, dc = d, xc(n);

    // Forward elimination
    for (int i = 1; i < n; i++) {
        double mc = ac[i-1] / bc[i-1];
        bc[i] = bc[i] - mc * cc[i-1];
        dc[i] = dc[i] - mc * dc[i-1];
    }

    // Back substitution
    xc[n-1] = dc[n-1] / bc[n-1];
    for (int i = n-2; i >= 0; i--) {
        xc[i] = (dc[i] - cc[i] * xc[i+1]) / bc[i];
    }

    return xc;
}

int main() {
    int n = 100;

    // Initialize the main diagonal, super-diagonal, sub-diagonal
    vector<double> main_diag(n, 2.0004);
    vector<double> super_diag(n-1, -1);
    vector<double> sub_diag(n-1, -1);

    // Adjust the last sub-diagonal for boundary conditions
    sub_diag[n-2] = -2;

    // Initialize the RHS vector
    vector<double> rhs(n, 0.0);
    rhs[0] = 1;

    // Call Thomas algorithm function
    vector<double> result = TDMAsolver(sub_diag, main_diag, super_diag, rhs);

    // Output result
    for (double x : result) {
        cout << x << " ";
    }
    cout << endl;

    return 0;
}
