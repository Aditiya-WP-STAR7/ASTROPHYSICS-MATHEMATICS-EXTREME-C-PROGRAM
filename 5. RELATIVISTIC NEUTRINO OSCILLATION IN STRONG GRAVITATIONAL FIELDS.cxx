#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

// =======================
// Physical Constants
// =======================
const double PI = 3.141592653589793;
const double G  = 6.67430e-11;
const double c  = 3.0e8;

// =======================
// Complex Matrix (2x2)
// =======================
struct Matrix2 {
    complex<double> m[2][2];

    Matrix2() {
        m[0][0] = m[1][1] = 0.0;
        m[0][1] = m[1][0] = 0.0;
    }

    static Matrix2 identity() {
        Matrix2 I;
        I.m[0][0] = I.m[1][1] = 1.0;
        return I;
    }
};

// =======================
// Matrix Operations
// =======================
Matrix2 operator+(const Matrix2& A, const Matrix2& B) {
    Matrix2 R;
    for(int i=0;i<2;i++)
        for(int j=0;j<2;j++)
            R.m[i][j] = A.m[i][j] + B.m[i][j];
    return R;
}

Matrix2 operator-(const Matrix2& A, const Matrix2& B) {
    Matrix2 R;
    for(int i=0;i<2;i++)
        for(int j=0;j<2;j++)
            R.m[i][j] = A.m[i][j] - B.m[i][j];
    return R;
}

Matrix2 operator*(const Matrix2& A, const Matrix2& B) {
    Matrix2 R;
    for(int i=0;i<2;i++)
        for(int j=0;j<2;j++)
            for(int k=0;k<2;k++)
                R.m[i][j] += A.m[i][k] * B.m[k][j];
    return R;
}

Matrix2 operator*(complex<double> s, const Matrix2& A) {
    Matrix2 R;
    for(int i=0;i<2;i++)
        for(int j=0;j<2;j++)
            R.m[i][j] = s * A.m[i][j];
    return R;
}

// =======================
// Commutator [H,ρ]
// =======================
Matrix2 commutator(const Matrix2& H, const Matrix2& rho) {
    return H * rho - rho * H;
}

// =======================
// Gravitational Redshift
// =======================
double gravitationalRedshift(double M, double r) {
    return sqrt(1.0 - (2.0 * G * M) / (r * c * c));
}

// =======================
// Hamiltonian Construction
// =======================
Matrix2 hamiltonian(double dm2, double E, double theta, double gravFactor) {
    Matrix2 H;

    double delta = dm2 / (2.0 * E) * gravFactor;

    H.m[0][0] = -delta * cos(2.0 * theta);
    H.m[1][1] =  delta * cos(2.0 * theta);

    H.m[0][1] =  delta * sin(2.0 * theta);
    H.m[1][0] =  delta * sin(2.0 * theta);

    return H;
}

// =======================
// RK-Nyström Integration
// =======================
Matrix2 evolveDensityMatrix(
    const Matrix2& rho0,
    const Matrix2& H,
    double dLambda
) {
    complex<double> I(0.0, 1.0);

    Matrix2 k1 = (-I) * commutator(H, rho0);
    Matrix2 k2 = (-I) * commutator(H, rho0 + (0.5 * dLambda) * k1);
    Matrix2 k3 = (-I) * commutator(H, rho0 + (0.5 * dLambda) * k2);
    Matrix2 k4 = (-I) * commutator(H, rho0 + dLambda * k3);

    Matrix2 rho = rho0 +
        (dLambda / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);

    return rho;
}

// =======================
// Main CLI Program
// =======================
int main() {
    cout << fixed << setprecision(6);

    while(true) {
        cout << "\n=============================================\n";
        cout << " Relativistic Neutrino Oscillation Simulator\n";
        cout << " Density Matrix + Strong Gravitational Field\n";
        cout << "=============================================\n";

        double dm2, E, theta;
        double M, r;
        int steps;

        cout << "Mass-squared difference Δm² (eV²): ";
        cin >> dm2;

        cout << "Neutrino energy E (eV): ";
        cin >> E;

        cout << "Mixing angle θ (radians): ";
        cin >> theta;

        cout << "Central mass M (kg): ";
        cin >> M;

        cout << "Radial distance r (meters): ";
        cin >> r;

        cout << "Integration steps: ";
        cin >> steps;

        double grav = gravitationalRedshift(M, r);
        double dLambda = 1.0 / steps;

        Matrix2 rho;
        rho.m[0][0] = 1.0;
        rho.m[1][1] = 0.0;

        Matrix2 H = hamiltonian(dm2, E, theta, grav);

        for(int i=0;i<steps;i++) {
            rho = evolveDensityMatrix(rho, H, dLambda);
        }

        cout << "\nFinal Density Matrix ρ:\n";
        cout << "ρ11 = " << rho.m[0][0] << "\n";
        cout << "ρ12 = " << rho.m[0][1] << "\n";
        cout << "ρ21 = " << rho.m[1][0] << "\n";
        cout << "ρ22 = " << rho.m[1][1] << "\n";

        cout << "\nElectron Neutrino Survival Probability:\n";
        cout << abs(rho.m[0][0]) << "\n";

        char repeat;
        cout << "\nRun another calculation? (y/n): ";
        cin >> repeat;

        if(repeat != 'y' && repeat != 'Y')
            break;
    }

    cout << "\nSimulation terminated. Stay curious.\n";
    return 0;
}
