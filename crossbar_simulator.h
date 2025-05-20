#ifndef CROSSBAR_SIMULATOR_H_
#define CROSSBAR_SIMULATOR_H_

#include "memristor_model/memristor.h"
#include "crossbar_model/linear_crossbar_solver.h"

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

    CrossbarSimulator(int M, int N) : linear_solver() {
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

    template <typename MemristorType>
    void Initialize() {
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                RRAM[i][j] = std::make_unique<MemristorType>();
            }
        }
    }
    
    void SetRRAM(std::vector<std::vector<bool>> weights);
    Eigen::VectorXf NonlinearSolve(
        Eigen::VectorXf Vguess,
        const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
        const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
        std::string method = "fixed-point"
    );
    std::vector<std::vector<float>> ApplyVoltage(
        Eigen::VectorXf Vguess,
        const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
        const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
        float dt, std::string method = "fixed-point"
    );
    std::vector<float> CalculateIout(Eigen::VectorXf Vout);

    void SetCrossbarParameters(float Rswl1, float Rswl2, float Rsbl1, float Rsbl2, float Rwl, float Rbl);
    void UpdatePrecomputeG_ABCD();
};

#endif  // CROSSBAR_SIMULATOR_H_
