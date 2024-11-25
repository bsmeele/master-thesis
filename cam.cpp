#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <chrono>

std::vector<float> solve_cam(
    const Eigen::MatrixXf& G,
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
    const float Rswl1, const float Rswl2, const float Rsbl1, const float Rsbl2,
    const float Rwl, const float Rbl,
    const bool print = false
    ) {
    int M = G.rows();
    int N = G.cols();

    float Gswl1 = 1/Rswl1;
    float Gswl2 = 0;
    float Gsbl1 = 0;
    float Gsbl2 = 1/Rsbl2;

    float Gwl = 1/Rwl;
    float Gbl = 1/Rbl;

    Eigen::SparseMatrix<float> G_ABCD(2*M*N, 2*M*N);

    // Submatrix A
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            if (j == 0) {
                G_ABCD.insert(i*N + j, i*N + j) = Gswl1 + G(i, j) + Gwl;
                G_ABCD.insert(i*N + j, i*N + j+1) = -Gwl;
            } else if (j == N-1) {
                G_ABCD.insert(i*N + j, i*N + j) = Gswl2 + G(i, j) + Gwl;
                G_ABCD.insert(i*N + j, i*N + j-1) = -Gwl;
            } else {
                G_ABCD.insert(i*N + j, i*N+j) = G(i, j) + 2*Gwl;
                G_ABCD.insert(i*N + j, i*N + j-1) = -Gwl;
                G_ABCD.insert(i*N + j, i*N + j+1) = -Gwl;
            }
        }
    }

    // Submatrix B
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            G_ABCD.insert(i*N + j, i*N + j + M*N) = -G(i, j);
        }
    }

    // Submatrix C
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            G_ABCD.insert(i*N + j + M*N, i*N + j) = G(i, j);
        }
    }

    // Submatrix D
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            if (i == 0) {
                G_ABCD.insert(i*N + j + M*N, i*N + j + M*N) = -Gsbl1 + -G(i, j) + -Gbl;
                G_ABCD.insert(i*N + j + M*N, (i+1)*N + j + M*N) = Gbl;
            } else if (i == M-1) {
                G_ABCD.insert(i*N + j + M*N, i*N + j + M*N) = -Gsbl2 + -G(i, j) + -Gbl;
                G_ABCD.insert(i*N + j + M*N, (i-1)*N + j + M*N) = Gbl;
            } else {
                G_ABCD.insert(i*N + j + M*N, i*N + j + M*N) = -G(i, j) + -2*Gbl;
                G_ABCD.insert(i*N + j + M*N, (i-1)*N + j + M*N) = Gbl;
                G_ABCD.insert(i*N + j + M*N, (i+1)*N + j + M*N) = Gbl;
            }
        }
    }

    if (print) {
        std::cout << "G_ABCD:\n" << G_ABCD << std::endl << std::endl;
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
    
    if (print) {
        std::cout << "E:\n" << E << std::endl << std::endl;
    }

    // Eigen::SparseLU<Eigen::SparseMatrix<float>> solver;
    // Eigen::ConjugateGradient<Eigen::SparseMatrix<float>> solver;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<float>> solver;

    solver.compute(G_ABCD);
    Eigen::VectorXf V = solver.solve(E);

    if (print) {
        std::cout << "Va:\n" << V.head(M*N) << std::endl << std::endl;
        std::cout << "Vb:\n" << V.tail(M*N) << std::endl << std::endl;
        std::cout << "V:\n" << V.head(M*N) - V.tail(M*N) << std::endl << std::endl;
    }

    // Calculate Iout
    std::vector<float> Iout;
    for (int j = 0; j < N; j++) {
        float Ioutj = 0;
        for (int i = 0; i < M; i++) {
            Ioutj += (V(i*N + j) - V(i*N + j + M*N)) * G(i,j);
        }
        Iout.push_back(Ioutj);
    }

    return Iout;
}

int main(int argc, char* argv[]) {
    srand((unsigned int) time(0));

    int M = (argc >= 2) ? std::atoi(argv[1]) : 16;
    int N = (argc >= 3) ? std::atoi(argv[2]) : 16;

    int runs = (argc >= 4) ? std::atoi(argv[3]) : 1;

    bool print = false;
    if (std::atoi(argv[4]) == 1) { print = true; }

    long long total_time = 0;

    for (int i = 0; i < runs; i++) {
        float Rmin = 100.;
        float Rmax = 1000.;
        float Vdd = 5.;

        float Rswl1 = 3.;
        float Rswl2 = INFINITY;
        float Rsbl1 = INFINITY;
        float Rsbl2 = 5.;

        float Rwl = 3.;
        float Rbl = 2.;

        if (print) {
            std::cout << "Rswl1: " << Rswl1 << std::endl;
            // std::cout << "Rswl2: " << Rswl2 << std::endl;
            // std::cout << "Rsbl1: " << Rsbl1 << std::endl;
            std::cout << "Rsbl2: " << Rsbl2 << std::endl;

            std::cout << "Rwl: " << Rwl << std::endl;
            std::cout << "Rbl: " << Rbl << std::endl << std::endl;
        }

        Eigen::MatrixXf R = Eigen::MatrixXf::Random(M, N);
        R = (R + Eigen::MatrixXf::Constant(M, N, 1.0)) * 0.5 * (Rmax - Rmin) + Eigen::MatrixXf::Constant(M, N, Rmin);
        // R = Eigen::MatrixXf::Zero(M, N);
        R(0, 0) = 10.;
        R(0, 1) = 15.;
        R(0, 2) = 20.;
        R(1, 0) = 25.;
        R(1, 1) = 30.;
        R(1, 2) = 35.;
        R(2, 0) = 40.;
        R(2, 1) = 45.;
        R(2, 2) = 50.;

        if (print) {
            std::cout << "R:\n" << R << std::endl << std::endl;
        }

        Eigen::MatrixXf G = R.cwiseInverse();

        if (print) {
            std::cout << "G:\n" << G << std::endl << std::endl;
        }

        Eigen::VectorXf Vappwl1 = Eigen::VectorXf::Random(M);
        Vappwl1 = (Vappwl1.array() > 0.5).select(Eigen::VectorXf::Constant(M, Vdd), Eigen::VectorXf::Zero(M));
        Vappwl1(0) = 5;
        Vappwl1(1) = 7;
        Vappwl1(2) = 9;

        Eigen::VectorXf Vappwl2 = Eigen::VectorXf::Zero(M);
        Eigen::VectorXf Vappbl1 = Eigen::VectorXf::Zero(M);
        Eigen::VectorXf Vappbl2 = Eigen::VectorXf::Zero(M);

        if (print) {
            std::cout << "Vappwl1:\n" << Vappwl1 << std::endl << std::endl;
        }

        auto start_time = std::chrono::high_resolution_clock::now();

        std::vector<float> Iout = solve_cam(G, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, print);

        auto end_time = std::chrono::high_resolution_clock::now();

        if (print) {
            std::cout << "Iout:" << std::endl;
            for (int j = 0; j < N; j++) {
                std::cout << Iout[j] << std::endl;
            }
            std::cout << std::endl;
        }

        auto execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        std::cout << "Execution time: " << execution_time << " (ms)" << std::endl;

        total_time += execution_time;
    }
    
    std::cout << "Average execution time: " << total_time/runs << " ms" << std::endl;
}
