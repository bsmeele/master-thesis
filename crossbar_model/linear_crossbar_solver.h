#ifndef LINEAR_CROSSBAR_SOLVER_H_
#define LINEAR_CROSSBAR_SOLVER_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

// Computes the nodal voltages of a crossbar based on crossbar resistances and device resistances
Eigen::VectorXf SolveCam(
    const Eigen::MatrixXf& G,  // Conductance matrix of the devices in the crossbar
    const Eigen::VectorXf& Vguess,  // Initial guess for the nodal voltages. Supplying a zero vector acts as if no guess is given
    Eigen::SparseMatrix<float>& G_ABCD,  // A partially precomputed version of the G_ABCD matrix, obtaineed from the PartiallyPrecomputeG_ABCD() function
    const Eigen::VectorXf& E,  // The precomputed E vector, obtained from the ComputeE() function
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,  // Applied voltages to the wordlines of the crossbar
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,  // Applied voltages to the bitlines of the crossbar
    const float Rswl1, const float Rswl2, const float Rsbl1, const float Rsbl2,  // Resitances of the wordline and bitline voltage sources
    const float Rwl, const float Rbl,  // Wordline and bitline resistances
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float>>& solver,  // Solver that is used to solve the system
    const bool print = false  // Boolean variable to print some debug information, default false
);

// Partially precomputes some elements of the G_ABCD matrix used by the SolveCam() function
Eigen::SparseMatrix<float> PartiallyPrecomputeG_ABCD(
    const int M, const int N,  // Crossbar dimensions
    const float Rswl1, const float Rswl2, const float Rsbl1, const float Rsbl2,  // Crossbar access resistors
    const float Rwl, const float Rbl  // Wordline and bitline resistances
);

// Precomputes the E vector used by the SolveCam() function
Eigen::VectorXf ComputeE(
    const int M, const int N,  // Crossbar dimensions
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,  // Applied voltages to the wordlines of the crossbar
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,  // Applied voltages to the bitlines of the crossbar
    const float Rswl1, const float Rswl2, const float Rsbl1, const float Rsbl2,  // Resitances of the wordline and bitline voltage sources
    const float Rwl, const float Rbl  // Wordline and bitline resistances
);

#endif  // LINEAR_CROSSBAR_SOLVER_H_
