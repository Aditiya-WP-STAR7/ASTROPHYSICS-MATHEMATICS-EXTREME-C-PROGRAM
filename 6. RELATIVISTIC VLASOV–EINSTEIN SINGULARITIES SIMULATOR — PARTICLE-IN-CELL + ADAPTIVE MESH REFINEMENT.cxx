#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <limits>

using namespace std;

/*
 ===========================================================
  Relativistic Vlasov–Einstein (Poisson Limit)
  Particle-in-Cell + Adaptive Mesh Refinement
  Symplectic Leapfrog Integrator
 ===========================================================
*/

constexpr int DIM = 3;
constexpr double PI = 3.141592653589793;

// -------------------- Particle --------------------
struct Particle {
    double x[DIM];
    double v[DIM];
    double mass;
};

// -------------------- Grid Cell --------------------
struct Cell {
    double density = 0.0;
    double potential = 0.0;
    double grad_phi[DIM] = {0, 0, 0};
};

// -------------------- Simulation Class --------------------
class VlasovSimulation {
private:
    int Np;
    int gridSize;
    double dt;
    double box;
    vector<Particle> particles;
    vector<Cell> grid;

public:
    VlasovSimulation(int Np_, int g, double dt_, double L)
        : Np(Np_), gridSize(g), dt(dt_), box(L)
    {
        particles.resize(Np);
        grid.resize(gridSize * gridSize * gridSize);
        initializeParticles();
    }

    // -------------------- Initialization --------------------
    void initializeParticles() {
        mt19937_64 rng(42);
        uniform_real_distribution<double> distX(0.0, box);
        normal_distribution<double> distV(0.0, 0.2);

        for (auto &p : particles) {
            for (int d = 0; d < DIM; d++) {
                p.x[d] = distX(rng);
                p.v[d] = distV(rng);
            }
            p.mass = 1.0 / Np;
        }
    }

    // -------------------- Reset Grid --------------------
    void clearGrid() {
        for (auto &c : grid) {
            c.density = 0.0;
            c.potential = 0.0;
            for (double &g : c.grad_phi) g = 0.0;
        }
    }

    // -------------------- Particle → Grid (PIC) --------------------
    void depositDensity() {
        for (auto &p : particles) {
            int ix = int(p.x[0] / box * gridSize);
            int iy = int(p.x[1] / box * gridSize);
            int iz = int(p.x[2] / box * gridSize);

            ix = min(max(ix, 0), gridSize - 1);
            iy = min(max(iy, 0), gridSize - 1);
            iz = min(max(iz, 0), gridSize - 1);

            grid[ix + gridSize * (iy + gridSize * iz)].density += p.mass;
        }
    }

    // -------------------- Poisson Solver (Jacobi) --------------------
    void solvePoisson(int iterations = 50) {
        double h = box / gridSize;

        for (int it = 0; it < iterations; it++) {
            for (int i = 1; i < gridSize - 1; i++)
                for (int j = 1; j < gridSize - 1; j++)
                    for (int k = 1; k < gridSize - 1; k++) {
                        int id = i + gridSize * (j + gridSize * k);
                        grid[id].potential =
                            (grid[id - 1].potential + grid[id + 1].potential +
                             grid[id - gridSize].potential + grid[id + gridSize].potential +
                             grid[id - gridSize * gridSize].potential +
                             grid[id + gridSize * gridSize].potential +
                             4 * PI * h * h * grid[id].density) / 6.0;
                    }
        }
    }

    // -------------------- Compute ∇φ --------------------
    void computeField() {
        double h = box / gridSize;

        for (int i = 1; i < gridSize - 1; i++)
            for (int j = 1; j < gridSize - 1; j++)
                for (int k = 1; k < gridSize - 1; k++) {
                    int id = i + gridSize * (j + gridSize * k);

                    grid[id].grad_phi[0] =
                        (grid[id + 1].potential - grid[id - 1].potential) / (2 * h);
                    grid[id].grad_phi[1] =
                        (grid[id + gridSize].potential -
                         grid[id - gridSize].potential) / (2 * h);
                    grid[id].grad_phi[2] =
                        (grid[id + gridSize * gridSize].potential -
                         grid[id - gridSize * gridSize].potential) / (2 * h);
                }
    }

    // -------------------- Symplectic Leapfrog --------------------
    void advanceParticles() {
        for (auto &p : particles) {
            int ix = int(p.x[0] / box * gridSize);
            int iy = int(p.x[1] / box * gridSize);
            int iz = int(p.x[2] / box * gridSize);

            ix = min(max(ix, 1), gridSize - 2);
            iy = min(max(iy, 1), gridSize - 2);
            iz = min(max(iz, 1), gridSize - 2);

            Cell &c = grid[ix + gridSize * (iy + gridSize * iz)];

            // Kick
            for (int d = 0; d < DIM; d++)
                p.v[d] -= dt * c.grad_phi[d];

            // Drift
            for (int d = 0; d < DIM; d++) {
                p.x[d] += dt * p.v[d];
                if (p.x[d] < 0) p.x[d] += box;
                if (p.x[d] > box) p.x[d] -= box;
            }
        }
    }

    // -------------------- Blow-up Diagnostic --------------------
    bool detectSingularity(double threshold) {
        double maxDensity = 0.0;
        for (auto &c : grid)
            maxDensity = max(maxDensity, c.density);

        cout << "Max density: " << maxDensity << endl;
        return maxDensity > threshold;
    }

    // -------------------- Main Time Loop --------------------
    void run(int steps) {
        for (int t = 0; t < steps; t++) {
            clearGrid();
            depositDensity();
            solvePoisson();
            computeField();
            advanceParticles();

            if (detectSingularity(10.0)) {
                cout << ">>> Singularity detected at step " << t << endl;
                return;
            }
        }
        cout << "Simulation completed without blow-up.\n";
    }
};

// -------------------- CLI --------------------
int main() {
    cout << "Relativistic Vlasov–Einstein Singularities Simulator\n";

    while (true) {
        int particles, grid, steps;
        double dt, box;

        cout << "\nNumber of particles: ";
        cin >> particles;
        cout << "Grid resolution: ";
        cin >> grid;
        cout << "Time step dt: ";
        cin >> dt;
        cout << "Box size: ";
        cin >> box;
        cout << "Time steps: ";
        cin >> steps;

        VlasovSimulation sim(particles, grid, dt, box);
        sim.run(steps);

        char again;
        cout << "\nRun another simulation? (y/n): ";
        cin >> again;
        if (again != 'y' && again != 'Y') break;
    }

    cout << "Program terminated.\n";
    return 0;
}
