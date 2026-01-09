#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>

using namespace std;

using Complex = complex<double>;
const double PI = acos(-1.0);

/* ================= FFT IMPLEMENTATION ================= */

void fft(vector<Complex> &a, bool invert) {
    int n = a.size();
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j |= bit;
        if (i < j) swap(a[i], a[j]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * PI / len * (invert ? -1 : 1);
        Complex wlen(cos(ang), sin(ang));
        for (int i = 0; i < n; i += len) {
            Complex w(1);
            for (int j = 0; j < len / 2; j++) {
                Complex u = a[i + j];
                Complex v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }

    if (invert) {
        for (Complex &x : a)
            x /= n;
    }
}

/* ================= SIGNAL GENERATION ================= */

vector<double> generateInspiral(int N, double f0, double chirp) {
    vector<double> h(N);
    for (int i = 0; i < N; i++) {
        double t = i / (double)N;
        double phase = 2 * PI * (f0 * t + chirp * t * t);
        h[i] = sin(phase);
    }
    return h;
}

/* ================= INNER PRODUCT ================= */

double innerProduct(const vector<double> &a,
                    const vector<double> &b) {
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); i++)
        sum += a[i] * b[i];
    return sum;
}

/* ================= MATCHED FILTER ================= */

double matchedFilterSNR(const vector<double> &signal,
                        const vector<double> &templ) {

    int n = 1;
    while (n < (int)signal.size()) n <<= 1;

    vector<Complex> S(n), H(n);

    for (size_t i = 0; i < signal.size(); i++) {
        S[i] = signal[i];
        H[i] = templ[i];
    }

    fft(S, false);
    fft(H, false);

    for (int i = 0; i < n; i++)
        S[i] *= conj(H[i]);

    fft(S, true);

    vector<double> corr(signal.size());
    for (size_t i = 0; i < corr.size(); i++)
        corr[i] = abs(S[i]);

    double num = *max_element(corr.begin(), corr.end());
    double denom = sqrt(innerProduct(templ, templ));

    return num / denom;
}

/* ================= MAIN PROGRAM ================= */

int main() {
    cout << "=============================================\n";
    cout << " GRAVITATIONAL WAVE DETECTION - MATCHED FILTER\n";
    cout << " Binary Black Hole Inspiral Analysis (C++)\n";
    cout << "=============================================\n";

    while (true) {
        int N;
        double f0, chirp;

        cout << "\nEnter number of samples (power of 2 recommended): ";
        cin >> N;

        cout << "Enter initial frequency f0 (Hz): ";
        cin >> f0;

        cout << "Enter chirp mass factor: ";
        cin >> chirp;

        vector<double> signal = generateInspiral(N, f0, chirp);

        // Simulated noisy detector data
        vector<double> detector(N);
        for (int i = 0; i < N; i++)
            detector[i] = signal[i] + 0.3 * sin(20 * PI * i / N);

        double snr = matchedFilterSNR(detector, signal);

        cout << "\n>>> Optimal Matched Filter SNR = " << snr << "\n";

        cout << "\nRun another computation? (y/n): ";
        char choice;
        cin >> choice;
        if (choice != 'y' && choice != 'Y') break;
    }

    cout << "\nProgram finished. Physics never ends.\n";
    return 0;
}
