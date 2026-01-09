#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <iomanip>

using namespace std;

// =================== Physical Constants ===================
const double PI = 3.141592653589793;
const double EPSILON = 0.01;      // Slow-roll parameter
const double H = 1.0;             // Hubble scale (natural units)
const double LAMBDA = 0.05;       // Backreaction strength

// =================== Background Functions ===================
double scale_factor(double tau) {
    return -1.0 / (H * tau);
}

double z(double tau) {
    return scale_factor(tau) * sqrt(2.0 * EPSILON);
}

double zpp_over_z(double tau) {
    return 2.0 / (tau * tau);
}

// =================== Backreaction Source ===================
complex<double> backreaction(const complex<double>& v) {
    return LAMBDA * norm(v) * v;
}

// =================== Mukhanov-Sasaki Integrator ===================
void evolve_mode(
    double k,
    double tau_start,
    double tau_end,
    int steps,
    complex<double>& v,
    complex<double>& vp
) {
    double dt = (tau_end - tau_start) / steps;
    double tau = tau_start;

    for (int i = 0; i < steps; i++) {
        complex<double> vpp =
            - (k * k - zpp_over_z(tau)) * v
            + backreaction(v);

        // Semi-implicit stiff update
        vp += dt * vpp;
        v  += dt * vp;

        tau += dt;
    }
}

// =================== Power Spectrum ===================
double power_spectrum(double k, complex<double> v, double tau) {
    double zz = z(tau);
    double Rk = abs(v / zz);
    return (k * k * k) / (2.0 * PI * PI) * Rk * Rk;
}

// =================== MAIN PROGRAM ===================
int main() {
    cout << "\n=============================================\n";
    cout << " Cosmic Inflation Power Spectrum Simulator\n";
    cout << " Mukhanov–Sasaki with Nonlinear Backreaction\n";
    cout << "=============================================\n";

    while (true) {
        double k;
        cout << "\nEnter wavenumber k (e.g. 0.05): ";
        cin >> k;

        // Time domain
        double tau_i = -100.0 / k;   // Deep sub-horizon
        double tau_f = -0.001;       // Super-horizon
        int steps = 50000;

        // Bunch–Davies Initial Conditions
        complex<double> v =
            exp(complex<double>(0, -k * tau_i)) / sqrt(2.0 * k);
        complex<double> vp =
            complex<double>(0, -k) * v;

        evolve_mode(k, tau_i, tau_f, steps, v, vp);

        double Pk = power_spectrum(k, v, tau_f);

        cout << fixed << setprecision(8);
        cout << "\n---------------------------------------------\n";
        cout << " Power Spectrum Result\n";
        cout << " k           = " << k << "\n";
        cout << " P_R(k)      = " << Pk << "\n";
        cout << " ln(P_R)     = " << log(Pk) << "\n";
        cout << "---------------------------------------------\n";

        char repeat;
        cout << "\nRun another computation? (y/n): ";
        cin >> repeat;

        if (repeat != 'y' && repeat != 'Y') {
            cout << "\nSimulation terminated.\n";
            break;
        }
    }

    return 0;
}
