#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

Eigen::VectorXf solve_cam(
    const Eigen::MatrixXf& G,
    const Eigen::VectorXf& V_guess, Eigen::SparseMatrix<float> G_ABCD,
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
    const float Rswl1, const float Rswl2, const float Rsbl1, const float Rsbl2,
    const float Rwl, const float Rbl,
    const bool print = false
);

Eigen::SparseMatrix<float> partially_precompute_G_ABCD(
    const int M, const int N,
    const float Rswl1, const float Rswl2, const float Rsbl1, const float Rsbl2,
    const float Rwl, const float Rbl
);
