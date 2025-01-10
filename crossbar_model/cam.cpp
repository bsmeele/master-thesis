#include "cam.h"

#include <iostream>

// Based on:
// https://ieeexplore-ieee-org.tudelft.idm.oclc.org/document/6473873

Eigen::VectorXf solve_cam(
    const Eigen::MatrixXf& G,
    const Eigen::VectorXf& V_guess, Eigen::SparseMatrix<float>& G_ABCD,
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
    const float Rswl1, const float Rswl2, const float Rsbl1, const float Rsbl2,
    const float Rwl, const float Rbl,
    const bool print
    ) {
    // auto start_time = std::chrono::high_resolution_clock::now();

    int M = G.rows();
    int N = G.cols();

    float Gswl1 = 1/Rswl1;
    float Gswl2 = 0;
    float Gsbl1 = 0;
    float Gsbl2 = 1/Rsbl2;

    float Gwl = 1/Rwl;
    float Gbl = 1/Rbl;

    // Eigen::SparseMatrix<float> G_ABCD(2*M*N, 2*M*N);

    // auto timestamp1 = std::chrono::high_resolution_clock::now();
    // auto execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp1 - start_time).count();
    // std::cout << "Setup time: " << execution_time << " (ms)" << std::endl;
    // timestamp1 = std::chrono::high_resolution_clock::now();

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

    // auto timestamp2 = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp2 - timestamp1).count();
    // std::cout << "Submatrix A time: " << execution_time << " (ms)" << std::endl;
    // timestamp2 = std::chrono::high_resolution_clock::now();

    // Submatrix B
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            // G_ABCD.insert(i*N + j, i*N + j + M*N) = -G(i, j);
            G_ABCD.coeffRef(i*N + j, i*N + j + M*N) += -G(i, j);
        }
    }

    // auto timestamp3 = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp3 - timestamp2).count();
    // std::cout << "Submatrix B time: " << execution_time << " (ms)" << std::endl;
    // timestamp3 = std::chrono::high_resolution_clock::now();

    // Submatrix C
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            // G_ABCD.insert(i*N + j + M*N, i*N + j) = G(i, j);
            G_ABCD.coeffRef(i*N + j + M*N, i*N + j) += G(i, j);
        }
    }

    // auto timestamp4 = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp4 - timestamp3).count();
    // std::cout << "Submatrix C time: " << execution_time << " (ms)" << std::endl;
    // timestamp4 = std::chrono::high_resolution_clock::now();

    // Submatrix D
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            G_ABCD.coeffRef(i*N + j + M*N, i*N + j + M*N) += -G(i, j);
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

    // auto timestamp5 = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp5 - timestamp4).count();
    // std::cout << "Submatrix D time: " << execution_time << " (ms)" << std::endl;
    // timestamp5 = std::chrono::high_resolution_clock::now();

    if (print) {
        std::cout << "G_ABCD:\n" << G_ABCD.toDense() << std::endl << std::endl;
    }

    // Make E
    Eigen::SparseVector<float> E(2*M*N);
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            if (j == 0) {
                E.insert(i*N + j) = Vappwl1(i) * Gswl1;
            } else if (j == N-1) {
                E.insert(i*N + j) = Vappwl2[i] * Gswl2;
            }
        }
    }
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            if (i == 0) {
                E.insert(i*N + j + M*N) = -Vappbl1(j) * Gsbl1;
            } else if (i == M-1) {
                E.insert(i*N + j + M*N) = -Vappbl2(j) * Gsbl2;
            }
        }
    }

    // auto timestamp6 = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp6 - timestamp5).count();
    // std::cout << "Matrix E time: " << execution_time << " (ms)" << std::endl;
    // timestamp6 = std::chrono::high_resolution_clock::now();
    
    if (print) {
        std::cout << "E:\n" << E.toDense() << std::endl << std::endl;
    }

    // Eigen::SparseLU<Eigen::SparseMatrix<float>> solver;
    // Eigen::ConjugateGradient<Eigen::SparseMatrix<float>> solver;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<float>> solver;

    solver.compute(G_ABCD);
    // Eigen::VectorXf V = solver.solve(E);
    Eigen::VectorXf V = solver.solveWithGuess(E.toDense(), V_guess);

    // auto timestamp7 = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp7 - timestamp6).count();
    // std::cout << "Solve time: " << execution_time << " (ms)" << std::endl;
    // timestamp7 = std::chrono::high_resolution_clock::now();

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

    // // auto end_time = std::chrono::high_resolution_clock::now();
    // // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - timestamp7).count();
    // // std::cout << "Iout time: " << execution_time << " (ms)" << std::endl;

    // return Iout;

    return V;
}

Eigen::SparseMatrix<float> partially_precompute_G_ABCD(
    const int M, const int N,
    const float Rswl1, const float Rswl2, const float Rsbl1, const float Rsbl2,
    const float Rwl, const float Rbl
) {
    float Gswl1 = 1/Rswl1;
    float Gswl2 = 0;
    float Gsbl1 = 0;
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
            G_ABCD.insert(i*N + j, i*N + j + M*N) = 0;
        }
    }
    // Submatrix C
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            G_ABCD.insert(i*N + j + M*N, i*N + j) = 0;
        }
    }
    // Submatrix D
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            if (i == 0) {
                G_ABCD.insert(i*N + j + M*N, i*N + j + M*N) = -Gsbl1 + -Gbl;
                G_ABCD.insert(i*N + j + M*N, (i+1)*N + j + M*N) = Gbl;
            } else if (i == M-1) {
                G_ABCD.insert(i*N + j + M*N, i*N + j + M*N) = -Gsbl2 + -Gbl;
                G_ABCD.insert(i*N + j + M*N, (i-1)*N + j + M*N) = Gbl;
            } else {
                G_ABCD.insert(i*N + j + M*N, i*N + j + M*N) = -2*Gbl;
                G_ABCD.insert(i*N + j + M*N, (i-1)*N + j + M*N) = Gbl;
                G_ABCD.insert(i*N + j + M*N, (i+1)*N + j + M*N) = Gbl;
            }
        }
    }

    return G_ABCD;
}
