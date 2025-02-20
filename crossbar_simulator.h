#ifndef CROSSBAR_SIMULATOR_H_
#define CROSSBAR_SIMULATOR_H_

#include "memristor_model/JART_VCM_v1b_var.h"
#include "crossbar_model/linear_crossbar_solver.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <string>

class CrossbarSimulator {
    public:
    int M;
    int N;

    std::vector<std::vector<JART_VCM_v1b_var>> RRAM;
    std::vector<std::vector<bool>> access_transistors;

    float Rswl1;
    float Rswl2;
    float Rsbl1;
    float Rsbl2;
    float Rwl;
    float Rbl;

    Eigen::SparseMatrix<float> partial_G_ABCD;

    CrossbarSimulator(int M, int N) {
        M = M;
        N = N;

        // Initialize RRAM
        RRAM = std::vector<std::vector<JART_VCM_v1b_var>>(M, std::vector<JART_VCM_v1b_var>(N, JART_VCM_v1b_var()));

        // Initialize access transistors to all on
        access_transistors =  std::vector<std::vector<bool>>(M, std::vector<bool>(N, true));

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
    
    void SetRRAM(std::vector<std::vector<bool>> weights);
    // Eigen::VectorXf LinearSolve(
    //     Eigen::VectorXf Vguess,
    //     const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
    //     const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2
    // );  // Would need some assumptions to be able to make the G matrix
    std::vector<float> NonlinearSolve(
        Eigen::VectorXf Vguess,
        const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
        const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
        std::string method = "fixed-point"
    );
    std::vector<float> ApplyVoltage(
        Eigen::VectorXf Vguess,
        const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
        const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
        float dt, std::string method = "fixed-point"
    );
};

#endif  // CROSSBAR_SIMULATOR_H_
