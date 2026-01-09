#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>

using namespace std;

/* ============================================
   Fundamental Constants (Geometrized Units)
============================================ */
constexpr double c = 1.0;
constexpr double G = 1.0;

/* ============================================
   Kerr Black Hole Parameters
============================================ */
struct KerrMetric {
    double M;  // Black hole mass
    double a;  // Spin parameter
};

/* ============================================
   Fluid State in Relativistic Hydrodynamics
============================================ */
struct FluidState {
    double rho;        // Rest-mass density
    double pressure;   // Pressure
    double epsilon;    // Specific internal energy
    vector<double> u;  // 4-velocity (u^t, u^r, u^θ, u^φ)
};

/* ============================================
   Boyer–Lindquist Metric Tensor
============================================ */
vector<vector<double>> metricTensor(double r, double theta, const KerrMetric& k) {
    double Sigma = r*r + k.a*k.a * cos(theta)*cos(theta);
    double Delta = r*r - 2*k.M*r + k.a*k.a;

    vector<vector<double>> g(4, vector<double>(4, 0.0));

    g[0][0] = -(1.0 - 2.0*k.M*r/Sigma);
    g[1][1] = Sigma / Delta;
    g[2][2] = Sigma;
    g[3][3] = sin(theta)*sin(theta) *
              (r*r + k.a*k.a + 2*k.M*r*k.a*k.a*sin(theta)*sin(theta)/Sigma);

    g[0][3] = g[3][0] = -2*k.M*r*k.a*sin(theta)*sin(theta)/Sigma;

    return g;
}

/* ============================================
   Relativistic Enthalpy
============================================ */
double enthalpy(const FluidState& f) {
    return 1.0 + f.epsilon + f.pressure / f.rho;
}

/* ============================================
   Alpha-Viscosity Model
============================================ */
double alphaViscosity(double alpha, double soundSpeed, double scaleHeight) {
    return alpha * soundSpeed * scaleHeight;
}

/* ============================================
   Relativistic Shear Tensor (Simplified DNS)
============================================ */
vector<vector<double>> shearTensor(const FluidState& f) {
    vector<vector<double>> sigma(4, vector<double>(4, 0.0));

    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            sigma[i][j] = 0.5 * (f.u[i] * f.u[j]);

    return sigma;
}

/* ============================================
   Stress–Energy Tensor
============================================ */
vector<vector<double>> stressEnergyTensor(
        const FluidState& f,
        const vector<vector<double>>& g,
        double nu) {

    vector<vector<double>> T(4, vector<double>(4, 0.0));
    auto sigma = shearTensor(f);
    double h = enthalpy(f);

    for (int mu = 0; mu < 4; ++mu) {
        for (int nu_i = 0; nu_i < 4; ++nu_i) {
            T[mu][nu_i] =
                f.rho * h * f.u[mu] * f.u[nu_i]
                + f.pressure * g[mu][nu_i]
                - 2.0 * nu * sigma[mu][nu_i];
        }
    }
    return T;
}

/* ============================================
   Spectral Element Step (Conceptual DNS Core)
============================================ */
void spectralElementStep(FluidState& f, double dt) {
    for (int i = 0; i < 4; ++i) {
        f.u[i] += -dt * 0.01 * f.u[i]; // Dissipative turbulent decay
    }
}

/* ============================================
   Main Simulation Loop
============================================ */
int main() {
    cout << fixed << setprecision(6);

    while (true) {
        cout << "\n=== Relativistic Turbulent Accretion Disk Simulator ===\n";

        KerrMetric k;
        FluidState f;

        cout << "Black hole mass M: ";
        cin >> k.M;
        cout << "Spin parameter a (0 ≤ a ≤ M): ";
        cin >> k.a;

        cout << "Density rho: ";
        cin >> f.rho;
        cout << "Pressure p: ";
        cin >> f.pressure;
        cout << "Internal energy epsilon: ";
        cin >> f.epsilon;

        f.u = {1.0, 0.01, 0.0, 0.3};

        double r, theta;
        cout << "Radial coordinate r: ";
        cin >> r;
        cout << "Polar angle theta (rad): ";
        cin >> theta;

        double alpha, cs, H;
        cout << "Alpha viscosity parameter α: ";
        cin >> alpha;
        cout << "Sound speed c_s: ";
        cin >> cs;
        cout << "Disk scale height H: ";
        cin >> H;

        double nu = alphaViscosity(alpha, cs, H);
        auto g = metricTensor(r, theta, k);

        double dt = 0.01;
        int steps = 100;

        for (int i = 0; i < steps; ++i) {
            spectralElementStep(f, dt);
        }

        auto T = stressEnergyTensor(f, g, nu);

        cout << "\nStress–Energy Tensor T^{μν}:\n";
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                cout << setw(14) << T[i][j] << " ";
            }
            cout << "\n";
        }

        cout << "\nRun another simulation? (y/n): ";
        char choice;
        cin >> choice;
        if (choice != 'y' && choice != 'Y') break;
    }

    cout << "\nSimulation finished. Results ready for analysis.\n";
    return 0;
}
