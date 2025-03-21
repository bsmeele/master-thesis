#ifndef LINEAR_CROSSBAR_SOLVER_H_
#define LINEAR_CROSSBAR_SOLVER_H_

#include </home/earapidis/dev/eigen/Eigen/Dense>
#include </home/earapidis/dev/eigen/Eigen/Sparse>
#include <vector>

Eigen::VectorXf SolveCam(
    const Eigen::MatrixXf& G,
    const Eigen::VectorXf& Vguess, Eigen::SparseMatrix<float> G_ABCD,
    const Eigen::VectorXf& E,
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
    const float Rswl1, const float Rswl2, const float Rsbl1, const float Rsbl2,
    const float Rwl, const float Rbl,
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float>>& solver,
    const bool print = false
);

Eigen::SparseMatrix<float> PartiallyPrecomputeG_ABCD(
    const int M, const int N,
    const float Rswl1, const float Rswl2, const float Rsbl1, const float Rsbl2,
    const float Rwl, const float Rbl
);

Eigen::VectorXf ComputeE(
    const int M, const int N,
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,  // Applied voltages to the wordlines of the crossbar
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,  // Applied voltages to the bitlines of the crossbar
    const float Rswl1, const float Rswl2, const float Rsbl1, const float Rsbl2,
    const float Rwl, const float Rbl
);

#endif  // LINEAR_CROSSBAR_SOLVER_H_
