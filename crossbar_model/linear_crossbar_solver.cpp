#include "linear_crossbar_solver.h"

#include <iostream>

// Based on:
// https://ieeexplore-ieee-org.tudelft.idm.oclc.org/document/6473873

// Computes the nodal voltages of a crossbar based on crossbar resistances and device resistances
// For a mathemtical explanation of this method, see the reference paper
Eigen::VectorXf SolveCam(
    const Eigen::MatrixXf& G,  // Conductance matrix of the devices in the crossbar
    const Eigen::VectorXf& Vguess,  // Initial guess for the nodal voltages. Supplying a zero vector acts as if no guess is given
    Eigen::SparseMatrix<float> G_ABCD,  // A partially precomputed version of the G_ABCD matrix. This precompution is done with the PartiallyPrecomputeG_ABCD() function
    const Eigen::VectorXf& E,
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,  // Applied voltages to the wordlines of the crossbar
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,  // Applied voltages to the bitlines of the crossbar
    const float Rswl1, const float Rswl2, const float Rsbl1, const float Rsbl2,  // Resitances of the wordline and bitline voltage sources
    const float Rwl, const float Rbl,  // Wordline and bitline resistances of the crossbar
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float>>& solver,
    const bool print  // Boolean variable to print some debug information, default false
    ) {

    int M = G.rows();
    int N = G.cols();

    float Gswl1 = 1/Rswl1;
    float Gswl2 = 1/Rswl2;
    float Gsbl1 = 1/Rsbl1;
    float Gsbl2 = 1/Rsbl2;

    float Gwl = 1/Rwl;
    float Gbl = 1/Rbl;

    // Submatrix A
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            G_ABCD.coeffRef(i*N + j, i*N + j) += G(i, j);
            // if (j == 0) {
            //     G_ABCD.insert(i*N + j, i*N + j) = Gswl1 + G(i, j) + Gwl;
            //     G_ABCD.insert(i*N + j, i*N + j+1) = -Gwl;
            // } else if (j == N-1) {
            //     G_ABCD.insert(i*N + j, i*N + j) = Gswl2 + G(i, j) + Gwl;
            //     G_ABCD.insert(i*N + j, i*N + j-1) = -Gwl;
            // } else {
            //     G_ABCD.insert(i*N + j, i*N+j) = G(i, j) + 2*Gwl;
            //     G_ABCD.insert(i*N + j, i*N + j-1) = -Gwl;
            //     G_ABCD.insert(i*N + j, i*N + j+1) = -Gwl;
            // }
        }
    }

    // Submatrix B
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            // G_ABCD.insert(i*N + j, i*N + j + M*N) = -G(i, j);
            G_ABCD.coeffRef(i*N + j, i*N + j + M*N) = -G(i, j);
        }
    }

    // Submatrix C
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            // G_ABCD.insert(i*N + j + M*N, i*N + j) = G(i, j);
            // G_ABCD.coeffRef(i*N + j + M*N, i*N + j) += G(i, j);
            G_ABCD.coeffRef(i*N + j + M*N, i*N + j) = -G(i, j);
        }
    }

    // Submatrix D
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            // G_ABCD.coeffRef(i*N + j + M*N, i*N + j + M*N) += -G(i, j);
            G_ABCD.coeffRef(i*N + j + M*N, i*N + j + M*N) += G(i, j);
            // if (i == 0) {
            //     G_ABCD.insert(i*N + j + M*N, i*N + j + M*N) = -Gsbl1 + -G(i, j) + -Gbl;
            //     G_ABCD.insert(i*N + j + M*N, (i+1)*N + j + M*N) = Gbl;
            // } else if (i == M-1) {
            //     G_ABCD.insert(i*N + j + M*N, i*N + j + M*N) = -Gsbl2 + -G(i, j) + -Gbl;
            //     G_ABCD.insert(i*N + j + M*N, (i-1)*N + j + M*N) = Gbl;
            // } else {
            //     G_ABCD.insert(i*N + j + M*N, i*N + j + M*N) = -G(i, j) + -2*Gbl;
            //     G_ABCD.insert(i*N + j + M*N, (i-1)*N + j + M*N) = Gbl;
            //     G_ABCD.insert(i*N + j + M*N, (i+1)*N + j + M*N) = Gbl;
            // }
        }
    }

    if (print) {
        std::cout << "G_ABCD:\n" << G_ABCD.toDense() << std::endl << std::endl;
    }
    
    if (print) {
        std::cout << "E:\n" << E << std::endl << std::endl;
    }

    // Eigen::SparseLU<Eigen::SparseMatrix<float>> solver;
    // Eigen::ConjugateGradient<Eigen::SparseMatrix<float>> solver;
    // Eigen::BiCGSTAB<Eigen::SparseMatrix<float>> solver;

    solver.compute(G_ABCD);
    // Eigen::VectorXf V = solver.solve(E);
    Eigen::VectorXf V = solver.solveWithGuess(E, Vguess);

    if (print) {
        std::cout << "Va:\n" << V.head(M*N) << std::endl << std::endl;
        std::cout << "Vb:\n" << V.tail(M*N) << std::endl << std::endl;
        std::cout << "V:\n" << V.head(M*N) - V.tail(M*N) << std::endl << std::endl;
    }

    // Calculate Iout
    // std::vector<float> Iout;
    // for (int j = 0; j < N; j++) {
    //     float Ioutj = 0;
    //     for (int i = 0; i < M; i++) {
    //         Ioutj += (V(i*N + j) - V(i*N + j + M*N)) * G(i,j);
    //     }
    //     Iout.push_back(Ioutj);
    // }

    // return Iout;

    return V;
}

// Partially precomputes some elements of the G_ABCD matrix used by the SolveCam() function
Eigen::SparseMatrix<float> PartiallyPrecomputeG_ABCD(
    const int M, const int N,
    const float Rswl1, const float Rswl2, const float Rsbl1, const float Rsbl2,
    const float Rwl, const float Rbl
) {
    float Gswl1 = 1/Rswl1;
    float Gswl2 = 1/Rswl2;
    float Gsbl1 = 1/Rsbl1;
    float Gsbl2 = 1/Rsbl2;

    float Gwl = 1/Rwl;
    float Gbl = 1/Rbl;
        
    // Partially precompute G_ABCD based on the crossbar parasitics
    Eigen::SparseMatrix<float> G_ABCD(2*M*N, 2*M*N);
    // Submatrix A
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            if (j == 0) {
                G_ABCD.insert(i*N + j, i*N + j) = Gswl1 + Gwl;
                G_ABCD.insert(i*N + j, i*N + j+1) = -Gwl;
            } else if (j == N-1) {
                G_ABCD.insert(i*N + j, i*N + j) = Gswl2 + Gwl;
                G_ABCD.insert(i*N + j, i*N + j-1) = -Gwl;
            } else {
                G_ABCD.insert(i*N + j, i*N+j) = 2*Gwl;
                G_ABCD.insert(i*N + j, i*N + j-1) = -Gwl;
                G_ABCD.insert(i*N + j, i*N + j+1) = -Gwl;
            }
        }
    }
    // Submatrix B
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            G_ABCD.insert(i*N + j, i*N + j + M*N) = 1e-6;
        }
    }
    // Submatrix C
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            G_ABCD.insert(i*N + j + M*N, i*N + j) = 1e-6;
        }
    }
    // Submatrix D
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            if (i == 0) {
                // G_ABCD.insert(i*N + j + M*N, i*N + j + M*N) = -Gsbl1 + -Gbl;
                // G_ABCD.insert(i*N + j + M*N, (i+1)*N + j + M*N) = Gbl;
                G_ABCD.insert(i*N + j + M*N, i*N + j + M*N) = Gsbl1 + Gbl;
                G_ABCD.insert(i*N + j + M*N, (i+1)*N + j + M*N) = -Gbl;
            } else if (i == M-1) {
                // G_ABCD.insert(i*N + j + M*N, i*N + j + M*N) = -Gsbl2 + -Gbl;
                // G_ABCD.insert(i*N + j + M*N, (i-1)*N + j + M*N) = Gbl;
                G_ABCD.insert(i*N + j + M*N, i*N + j + M*N) = Gsbl2 + Gbl;
                G_ABCD.insert(i*N + j + M*N, (i-1)*N + j + M*N) = -Gbl;
            } else {
                // G_ABCD.insert(i*N + j + M*N, i*N + j + M*N) = -2*Gbl;
                // G_ABCD.insert(i*N + j + M*N, (i-1)*N + j + M*N) = Gbl;
                // G_ABCD.insert(i*N + j + M*N, (i+1)*N + j + M*N) = Gbl;
                G_ABCD.insert(i*N + j + M*N, i*N + j + M*N) = 2*Gbl;
                G_ABCD.insert(i*N + j + M*N, (i-1)*N + j + M*N) = -Gbl;
                G_ABCD.insert(i*N + j + M*N, (i+1)*N + j + M*N) = -Gbl;
            }
        }
    }

    return G_ABCD;
}

Eigen::VectorXf ComputeE(
    const int M, const int N,
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,  // Applied voltages to the wordlines of the crossbar
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,  // Applied voltages to the bitlines of the crossbar
    const float Rswl1, const float Rswl2, const float Rsbl1, const float Rsbl2,
    const float Rwl, const float Rbl
) {
    float Gswl1 = 1/Rswl1;
    float Gswl2 = 1/Rswl2;
    float Gsbl1 = 1/Rsbl1;
    float Gsbl2 = 1/Rsbl2;

    float Gwl = 1/Rwl;
    float Gbl = 1/Rbl;

    // Make E
    Eigen::VectorXf E = Eigen::VectorXf::Zero(2*M*N);
    // Ew
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            if (j == 0) {
                E(i*N + j) = Vappwl1(i) * Gswl1;
            } else if (j == N-1) {
                E(i*N + j) = Vappwl2[i] * Gswl2;
            }
        }
    }
    // Eb
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            if (i == 0) {
                // E.insert(i*N + j + M*N) = -Vappbl1(j) * Gsbl1;
                E(i*N + j + M*N) = Vappbl1(j) * Gsbl1;
            } else if (i == M-1) {
                // E.insert(i*N + j + M*N) = -Vappbl2(j) * Gsbl2;
                E(i*N + j + M*N) = Vappbl2(j) * Gsbl2;
            }
        }
    }

    return E;
}