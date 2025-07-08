#ifndef CROSSBAR_SIMULATOR_H_
#define CROSSBAR_SIMULATOR_H_

#include "memristor_model/memristor.h"
#include "crossbar_model/linear_crossbar_solver.h"
#include "threadpool.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <string>
#include <memory>

class CrossbarSimulator {
    public:
    int M;
    int N;

    std::vector<std::vector<std::unique_ptr<Memristor>>> RRAM;
    std::vector<std::vector<bool>> access_transistors;

    float Rswl1;
    float Rswl2;
    float Rsbl1;
    float Rsbl2;
    float Rwl;
    float Rbl;

    Eigen::SparseMatrix<float> partial_G_ABCD;

    Eigen::ConjugateGradient<Eigen::SparseMatrix<float>> linear_solver;

    ThreadPool pool;

    // Crossbar constructor. Should be followed by the Initialize() function to instantiate the correct memristor model
    CrossbarSimulator(int M, int N) : linear_solver(), pool(simulation_num_threads) {
        this->M = M;
        this->N = N;

        // Initialize RRAM
        // RRAM = std::vector<std::vector<std::unique_ptr<Memristor>>>(M, std::vector<std::unique_ptr<Memristor>>(N));
        RRAM.resize(M);
        for (int i = 0; i < M; ++i) {
            RRAM[i].resize(N);
        }
        // for (int i = 0; i < M; ++i) {
        //     for (int j = 0; j < N; ++j) {
        //         RRAM[i][j] = std::make_unique<MemristorType>();
        //     }
        // }

        // Initialize access transistors to all on
        access_transistors = std::vector<std::vector<bool>>(M, std::vector<bool>(N, true));

        // Initialize parasitics to some default
        Rswl1 = 1;
        Rswl2 = INFINITY;
        Rsbl1 = INFINITY;
        Rsbl2 = 1;
        Rwl = 1;
        Rbl = 1;

        // Precompute partial_G_ABCD
        partial_G_ABCD = PartiallyPrecomputeG_ABCD(M, N, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);
    }

    // Instantiates the crossbar with the provided memristor model
    template <typename MemristorType>
    void Initialize() {
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                RRAM[i][j] = std::make_unique<MemristorType>();
            }
        }
    }
    
    // Sets the state of every memristor in the crossbar
    // True sets to LRS, false sets to HRS
    void SetRRAM(std::vector<std::vector<bool>> weights);

    // Sets the access transistors of the memristors
    // Each element in the input vector sets all access transistors for one row
    void SetAccessTransistors(std::vector<bool> gate_lines);

    // Calculates the nodal voltages of the crossbar assuming the memristors are linear devices. Does not evolve memristor state
    Eigen::VectorXf LinearSolve(
        Eigen::VectorXf Vguess,
        const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
        const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2
    );

    // Calculates the nodal voltages of the crossbar assuming the memristors are non-linear devices. Does not evolve memristor state
    Eigen::VectorXf NonlinearSolve(
        Eigen::VectorXf Vguess,
        const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
        const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
        std::string method = "fixed-point"
    );

    // Simulates applying a voltage to the crossbar and returns the current running through each memristor. Includes evolving memristor state
    std::vector<std::vector<float>> ApplyVoltage(
        Eigen::VectorXf Vguess,
        const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
        const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
        const float dt,
        bool linear = false, std::string method = "fixed-point"
    );

    // Calculates the column current based on the nodal voltages of the crossbar
    std::vector<float> CalculateIout(const Eigen::VectorXf& Vout);

    float Energy(
        const Eigen::VectorXf& Vout,
        const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
        const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
        const float dt
    );

    // Encompases a complete simulation routine
    void Simulate(
        const std::vector<bool> Vwl1, const std::vector<bool> Vwl2,  // Applied voltages to the wordlines of the crossbar
        const std::vector<bool> Vbl1, const std::vector<bool> Vbl2,  // Applied voltages to the bitlines of the crossbar
        const std::vector<std::vector<bool>> weights,  // matrix of weights corresponding to each crossbar. Writing weights is not simulated, instead the SetRRAM() function is used
        const std::vector<std::array<float, 2>> waveform,  // Description of the waveform. Each element is a breakpoint consisting of a timestamp and a voltage. The wave is constructed by linearly interpolating between two breakpoints
        const float dt,  // Time step size used for simulation
        std::vector<std::vector<float>>& Iout,  // Output matrix for currents running through individual memristors. Will be cleared before use
        std::vector<float>& Iout_MAC,  // Output vector for column currents. Will be cleared before use
        bool linear = false, bool drift = true
    );

    // Sets the crossbar parameters and recalculates the precomputed G_ABCD matrix. Parameters updated through other means will not correctly be used in simulation
    void SetCrossbarParameters(float Rswl1, float Rswl2, float Rsbl1, float Rsbl2, float Rwl, float Rbl);


    void UpdatePrecomputeG_ABCD();
};

#endif  // CROSSBAR_SIMULATOR_H_
