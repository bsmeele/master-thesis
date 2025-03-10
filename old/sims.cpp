#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <chrono>

// Based on:
// https://ieeexplore-ieee-org.tudelft.idm.oclc.org/document/10559273

std::vector<float> solve_cam(
    const Eigen::MatrixXf& G,
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
    Eigen::VectorXf& V,
    const float Rswl1, const float Rswl2, const float Rsbl1, const float Rsbl2,
    const float Rwl, const float Rbl,
    const bool print = false
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

    // auto timestamp1 = std::chrono::high_resolution_clock::now();
    // auto execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp1 - start_time).count();
    // std::cout << "Setup time: " << execution_time << " (ms)" << std::endl;
    // timestamp1 = std::chrono::high_resolution_clock::now();

    // Submatrix A
    Eigen::MatrixXf A = Eigen::MatrixXf::Zero(M*N, M*N);
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            if (j == 0) {
                A(i*N + j, i*N + j) = Gswl1 + G(i, j) + Gwl;
                A(i*N + j, i*N + j+1) = -Gwl;
            } else if (j == N-1) {
                A(i*N + j, i*N + j) = Gswl2 + G(i, j) + Gwl;
                A(i*N + j, i*N + j-1) = -Gwl;
            } else {
                A(i*N + j, i*N+j) = G(i, j) + 2*Gwl;
                A(i*N + j, i*N + j-1) = -Gwl;
                A(i*N + j, i*N + j+1) = -Gwl;
            }
        }
    }

    if (print) {
        std::cout << "Submatrix A:\n" << A << std::endl << std::endl;
    }

    // auto timestamp2 = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp2 - timestamp1).count();
    // std::cout << "Submatrix A time: " << execution_time << " (ms)" << std::endl;
    // timestamp2 = std::chrono::high_resolution_clock::now();

    // Submatrix B
    Eigen::MatrixXf B = Eigen::MatrixXf::Zero(M*N, M*N);
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            B(i*N + j, i*N + j) = -G(i, j);
        }
    }

    if (print) {
        std::cout << "Submatrix B:\n" << B << std::endl << std::endl;
    }

    // auto timestamp3 = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp3 - timestamp2).count();
    // std::cout << "Submatrix B time: " << execution_time << " (ms)" << std::endl;
    // timestamp3 = std::chrono::high_resolution_clock::now();

    // Submatrix C
    Eigen::MatrixXf C = B;

    if (print) {
        std::cout << "Submatrix C:\n" << C << std::endl << std::endl;
    }

    // auto timestamp4 = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp4 - timestamp3).count();
    // std::cout << "Submatrix C time: " << execution_time << " (ms)" << std::endl;
    // timestamp4 = std::chrono::high_resolution_clock::now();

    // Submatrix D
    Eigen::MatrixXf D = Eigen::MatrixXf::Zero(M*N, M*N);
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            if (i == 0) {
                D(i*N + j, i*N + j) = Gsbl1 + G(i, j) + Gbl;
                D(i*N + j, (i+1)*N + j) = -Gbl;
            } else if (i == M-1) {
                D(i*N + j, i*N + j) = Gsbl2 + G(i, j) + Gbl;
                D(i*N + j, (i-1)*N + j) = -Gbl;
            } else {
                D(i*N + j, i*N + j) = G(i, j) + 2*Gbl;
                D(i*N + j, (i-1)*N + j) = -Gbl;
                D(i*N + j, (i+1)*N + j) = -Gbl;
            }
        }
    }

    if (print) {
        std::cout << "Submatrix D:\n" << D << std::endl << std::endl;
    }

    Eigen::MatrixXf G_ABCD(2*M*N, 2*M*N);
    G_ABCD.block(0, 0, A.rows(), A.cols()) = A;
    G_ABCD.block(0, M*N, B.rows(), B.cols()) = B;
    G_ABCD.block(M*N, 0, C.rows(), C.cols()) = C;
    G_ABCD.block(M*N, M*N, D.rows(), D.cols()) = D;

    // auto timestamp5 = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp5 - timestamp4).count();
    // std::cout << "Submatrix D time: " << execution_time << " (ms)" << std::endl;
    // timestamp5 = std::chrono::high_resolution_clock::now();

    if (print) {
        std::cout << "G_ABCD:\n" << G_ABCD << std::endl << std::endl;
    }

    // Make E
    Eigen::VectorXf E = Eigen::VectorXf::Zero(2*M*N);
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            if (j == 0) {
                E(i*N + j) = Vappwl1(i) * Gswl1;
            } else if (j == N-1) {
                E(i*N + j) = Vappwl2[i] * Gswl2;
            }
        }
    }
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            if (i == 0) {
                E(i*N + j + M*N) = -Vappbl1(j) * Gsbl1;
            } else if (i == M-1) {
                E(i*N + j + M*N) = -Vappbl2(j) * Gsbl2;
            }
        }
    }

    // auto timestamp6 = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp6 - timestamp5).count();
    // std::cout << "Matrix E time: " << execution_time << " (ms)" << std::endl;
    // timestamp6 = std::chrono::high_resolution_clock::now();
    
    if (print) {
        std::cout << "E:\n" << E << std::endl << std::endl;
    }

    // ----- SIMs decomposition -----
    float alpha = B.cwiseAbs().maxCoeff();

    // Calculate M
    Eigen::MatrixXf M_mat = Eigen::MatrixXf::Zero(2*M*N, 2*M*N);
    Eigen::MatrixXf M_top_left = A + (alpha * Eigen::MatrixXf::Identity(M*N, M*N) + B);
    Eigen::MatrixXf M_bottom_right = D + (alpha * Eigen::MatrixXf::Identity(M*N, M*N) + C);
    M_mat.block(0, 0, M_top_left.rows(), M_top_left.cols()) = M_top_left;
    M_mat.block(M*N, M*N, M_bottom_right.rows(), M_bottom_right.cols()) = M_bottom_right;

    // auto timestamp7 = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp7 - timestamp6).count();
    // std::cout << "M time: " << execution_time << " (ms)" << std::endl;
    // timestamp7 = std::chrono::high_resolution_clock::now();

    if (print) {
        std::cout << "M:\n" << M_mat << std::endl << std::endl;
    }

    Eigen::MatrixXf M_inv = M_mat.inverse();

    // auto timestamp8 = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp8 - timestamp7).count();
    // std::cout << "M inverse time: " << execution_time << " (ms)" << std::endl;
    // timestamp8 = std::chrono::high_resolution_clock::now();

    if (print) {
        std::cout << "M inverse:\n" << M_inv << std::endl << std::endl;
    }

    // Calculate N
    Eigen::MatrixXf N_mat(2*M*N, 2*M*N);
    Eigen::MatrixXf N_top_left = alpha * Eigen::MatrixXf::Identity(M*N, M*N) + B;
    Eigen::MatrixXf N_bottom_right = alpha * Eigen::MatrixXf::Identity(M*N, M*N) + C;
    N_mat.block(0, 0, N_top_left.rows(), N_top_left.cols()) = N_top_left;
    N_mat.block(0, M*N, B.rows(), B.cols()) = -B;
    N_mat.block(M*N, 0, C.rows(), C.cols()) = -C;
    N_mat.block(M*N, M*N, N_bottom_right.rows(), N_bottom_right.cols()) = N_bottom_right;

    // auto timestamp9 = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp9 - timestamp8).count();
    // std::cout << "N time: " << execution_time << " (ms)" << std::endl;
    // timestamp9 = std::chrono::high_resolution_clock::now();

    if (print) {
        std::cout << "N:\n" << N_mat << std::endl << std::endl;
    }

    // Calculate V
    V = M_inv * E;
    Eigen::MatrixXf prev_V;
    int max_it = 100;
    for (int i = 0; i < max_it; i++) {
        prev_V = V;
        V = M_inv * (N_mat * V + E);
        if ((V - prev_V).norm() < 1e-3) {
            std::cout << "Solved in " << i << " iterations" << std::endl << std::endl;
            break;
        } 
    }

    // auto timestamp10 = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp10 - timestamp9).count();
    // std::cout << "Solve time: " << execution_time << " (ms)" << std::endl;
    // timestamp10 = std::chrono::high_resolution_clock::now();

    // V = G_ABCD.colPivHouseholderQr().solve(E);

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

    // auto timestamp11 = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp11 - timestamp10).count();
    // std::cout << "Iout time: " << execution_time << " (ms)" << std::endl;

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
    Eigen::VectorXf V = Eigen::VectorXf::Zero(2*M*N);

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

        std::vector<float> Iout = solve_cam(G, Vappwl1, Vappwl2, Vappbl1, Vappbl2, V, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, print);

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
    
    if (runs > 1) {
        std::cout << "Average execution time: " << total_time/runs << " ms" << std::endl;
    }
}
