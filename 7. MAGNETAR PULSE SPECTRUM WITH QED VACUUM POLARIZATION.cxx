#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>

using namespace std;

/*
===========================================================
  MAGNETAR PULSE SPECTRUM WITH QED VACUUM POLARIZATION
  Euler–Heisenberg Effective Electrodynamics
  Author: (Your Name)
  Level : Research-Grade / MIT-Style Computational Physics
===========================================================
*/

// Physical constants (natural units where possible)
const double alpha = 1.0 / 137.035999;
const double me = 9.10938356e-31;          // electron mass (kg)
const double c = 2.99792458e8;             // speed of light (m/s)
const double hbar = 1.054571817e-34;
const double B_critical = 4.414e9;         // Schwinger critical field (Tesla)

// Numerical parameters
const int MAX_MODES = 1000;

// QED nonlinear correction coefficient
double qedCoefficient() {
    return (2.0 * alpha * alpha) / (45.0 * pow(me * c * c, 4));
}

// Effective refractive index due to vacuum polarization
double refractiveIndex(double B) {
    double ratio = B / B_critical;
    return 1.0 + qedCoefficient() * ratio * ratio;
}

// Compute eigenfrequencies of magnetar pulse modes
vector<double> computeSpectrum(double B, int modes) {
    vector<double> spectrum;
    double n_eff = refractiveIndex(B);

    for (int k = 1; k <= modes; ++k) {
        double omega = (k * M_PI * c) / n_eff;
        spectrum.push_back(omega);
    }
    return spectrum;
}

// Display spectrum
void displaySpectrum(const vector<double>& spectrum) {
    cout << "\n--- Magnetar Pulse Eigenmode Spectrum ---\n";
    cout << left << setw(10) << "Mode"
         << setw(25) << "Angular Frequency (rad/s)\n";
    cout << "----------------------------------------\n";

    for (size_t i = 0; i < spectrum.size(); ++i) {
        cout << setw(10) << i + 1
             << scientific << setw(25) << spectrum[i] << "\n";
    }
}

int main() {
    cout << "\n==============================================\n";
    cout << "  MAGNETAR PULSE SPECTRUM (QED VACUUM EFFECTS)\n";
    cout << "  Euler–Heisenberg Nonlinear Electrodynamics\n";
    cout << "==============================================\n";

    while (true) {
        double B;
        int modes;

        cout << "\nEnter magnetic field strength B (Tesla): ";
        cin >> B;

        cout << "Enter number of spectral modes (<= " << MAX_MODES << "): ";
        cin >> modes;

        if (cin.fail() || modes <= 0 || modes > MAX_MODES) {
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            cout << "Invalid input. Try again.\n";
            continue;
        }

        vector<double> spectrum = computeSpectrum(B, modes);
        displaySpectrum(spectrum);

        char choice;
        cout << "\nCompute another magnetar configuration? (y/n): ";
        cin >> choice;

        if (choice != 'y' && choice != 'Y') {
            cout << "\nExiting program. Continue pushing the frontier of physics.\n";
            break;
        }
    }

    return 0;
}
