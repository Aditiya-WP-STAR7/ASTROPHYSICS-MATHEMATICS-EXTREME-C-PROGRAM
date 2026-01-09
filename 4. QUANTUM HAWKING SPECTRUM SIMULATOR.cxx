#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

/*
============================================================
 Quantum Hawking Spectrum Simulator
 Lattice Quantum Gravity – Euclidean Path Integral
 Author: Research Prototype
============================================================
*/

constexpr double PI = 3.141592653589793;
constexpr double HBAR = 1.054571817e-34; // Reduced Planck constant
constexpr double G = 6.67430e-11;        // Gravitational constant
constexpr double C = 299792458.0;        // Speed of light
constexpr double KB = 1.380649e-23;      // Boltzmann constant

// Euclidean action for Schwarzschild black hole
double euclideanAction(double mass) {
    double beta = 8.0 * PI * G * mass / (C * C * C);
    double S = (beta * mass * C * C) / 2.0;
    return S;
}

// Quantum gravity correction expansion
double quantumCorrection(int order, double coefficient) {
    return coefficient * pow(HBAR, order);
}

// Lattice path integral approximation
double latticePathIntegral(
    double S_euclidean,
    const vector<double>& coefficients
) {
    double correctionSum = 0.0;

    for (size_t n = 0; n < coefficients.size(); ++n) {
        correctionSum += quantumCorrection(n + 1, coefficients[n]);
    }

    return exp(-(S_euclidean / HBAR) + correctionSum);
}

// Hawking temperature
double hawkingTemperature(double mass) {
    return (HBAR * C * C * C) / (8.0 * PI * G * mass * KB);
}

// Energy spectrum (Planck distribution)
double hawkingSpectrum(double energy, double temperature) {
    return energy / (exp(energy / (KB * temperature)) - 1.0);
}

int main() {
    cout << fixed << setprecision(6);

    while (true) {
        cout << "\n============================================\n";
        cout << " Quantum Hawking Radiation Spectrum Simulator\n";
        cout << " Lattice Quantum Gravity – CLI Research Tool\n";
        cout << "============================================\n";

        double mass;
        cout << "Enter Primordial Black Hole Mass (kg): ";
        cin >> mass;

        int latticeOrder;
        cout << "Enter Quantum Correction Order (N): ";
        cin >> latticeOrder;

        vector<double> coefficients(latticeOrder);
        for (int i = 0; i < latticeOrder; ++i) {
            cout << "Enter coefficient c_" << i + 1 << ": ";
            cin >> coefficients[i];
        }

        double S_E = euclideanAction(mass);
        double Z = latticePathIntegral(S_E, coefficients);
        double T = hawkingTemperature(mass);

        cout << "\n--- Results ---\n";
        cout << "Euclidean Action S_E       : " << S_E << "\n";
        cout << "Partition Function Z      : " << Z << "\n";
        cout << "Hawking Temperature (K)   : " << T << "\n";

        cout << "\nEnergy Spectrum Sample:\n";
        for (int i = 1; i <= 5; ++i) {
            double E = i * KB * T;
            cout << "E = " << E
                 << " -> Intensity = "
                 << hawkingSpectrum(E, T) << "\n";
        }

        cout << "\nRun another simulation? (y/n): ";
        char choice;
        cin >> choice;

        if (choice != 'y' && choice != 'Y') {
            cout << "\nExiting Quantum Gravity Simulator.\n";
            break;
        }
    }

    return 0;
}
