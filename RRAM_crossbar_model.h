#ifndef RRAM_CA_H_
#define RRAM_CA_H_

#include "memristor_model/JART_VCM_v1b_var.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vector>

Eigen::VectorXf BroydenInvSolve(
    std::vector<std::vector<JART_VCM_v1b_var>> RRAM,
    Eigen::VectorXf Vguess, Eigen::SparseMatrix<float> G_ABCD,
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
    const float Rswl1, const float Rswl2, const float Rsbl1, const float Rsbl2,
    const float Rwl, const float Rbl,
    const bool print = false
);

Eigen::VectorXf BroydenSolve(
    std::vector<std::vector<JART_VCM_v1b_var>> RRAM,
    Eigen::VectorXf Vguess, Eigen::SparseMatrix<float> G_ABCD,
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
    const float Rswl1, const float Rswl2, const float Rsbl1, const float Rsbl2,
    const float Rwl, const float Rbl,
    const bool print = false
);

Eigen::VectorXf NewtonRaphsonSolve(
    std::vector<std::vector<JART_VCM_v1b_var>> RRAM,
    Eigen::VectorXf Vguess, Eigen::SparseMatrix<float> G_ABCD,
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
    const float Rswl1, const float Rswl2, const float Rsbl1, const float Rsbl2,
    const float Rwl, const float Rbl,
    const bool print = false
);

Eigen::VectorXf FixedpointSolve(
    std::vector<std::vector<JART_VCM_v1b_var>> RRAM,
    Eigen::VectorXf Vguess, Eigen::SparseMatrix<float> G_ABCD,
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
    const float Rswl1, const float Rswl2, const float Rsbl1, const float Rsbl2,
    const float Rwl, const float Rbl,
    const bool print = false
);

#endif  // RRAM_CA_H_
