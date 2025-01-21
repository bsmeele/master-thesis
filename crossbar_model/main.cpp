#include "cam.h"

#include <iostream>
#include <chrono>
#include <cassert>

int main(int argc, char* argv[]) {
    srand((unsigned int) time(0));

    int M = (argc >= 2) ? std::atoi(argv[1]) : 16;
    int N = (argc >= 3) ? std::atoi(argv[2]) : 16;

    int runs = (argc >= 4) ? std::atoi(argv[3]) : 1;

    bool print = false;
    if (std::atoi(argv[4]) == 1) { print = true; }

    long long total_time = 0;
    Eigen::VectorXf V = Eigen::VectorXf::Zero(2*M*N);

    float Rswl1 = 3.;
    float Rswl2 = INFINITY;
    float Rsbl1 = INFINITY;
    float Rsbl2 = 5.;

    float Rwl = 3.;
    float Rbl = 2.;

    Eigen::SparseMatrix<float> G_ABCD = partially_precompute_G_ABCD(M, N, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);

    for (int i = 0; i < runs; i++) {
        float Rmin = 100.;
        float Rmax = 100000.;
        float Vdd = 1.5;

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
        // R(0, 0) = 10.;
        // R(0, 1) = 15.;
        // R(0, 2) = 20.;
        // R(1, 0) = 25.;
        // R(1, 1) = 30.;
        // R(1, 2) = 35.;
        // R(2, 0) = 40.;
        // R(2, 1) = 45.;
        // R(2, 2) = 50.;

        // R(0, 0) = 1e4;
        // R(0, 1) = 2e4;
        // R(0, 2) = 3e4;
        // R(1, 0) = 4e4;
        // R(1, 1) = 5e4;
        // R(1, 2) = 6e4;
        // R(2, 0) = 7e4;
        // R(2, 1) = 8e4;
        // R(2, 2) = 9e4;

        if (print) {
            std::cout << "R:\n" << R << std::endl << std::endl;
        }

        Eigen::MatrixXf G = R.cwiseInverse();

        if (print) {
            std::cout << "G:\n" << G << std::endl << std::endl;
        }

        Eigen::VectorXf Vappwl1 = Eigen::VectorXf::Random(M);
        Vappwl1 = (Vappwl1.array() > 0.5).select(Eigen::VectorXf::Constant(M, Vdd), Eigen::VectorXf::Zero(M));
        // Vappwl1(0) = 5;
        // Vappwl1(1) = 7;
        // Vappwl1(2) = 9;

        Eigen::VectorXf Vappwl2 = Eigen::VectorXf::Zero(M);
        Eigen::VectorXf Vappbl1 = Eigen::VectorXf::Zero(M);
        Eigen::VectorXf Vappbl2 = Eigen::VectorXf::Zero(M);

        if (print) {
            std::cout << "Vappwl1:\n" << Vappwl1 << std::endl << std::endl;
        }

        auto start_time = std::chrono::high_resolution_clock::now();

        Eigen::VectorXf Vout = solve_cam(G, V, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, false);

        auto end_time = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < Vout.size(); i++) {
            assert(!std::isnan(Vout(i)));
        }

        if (print) {
            std::vector<float> Iout;
            for (int j = 0; j < N; j++) {
                float Ioutj = 0;
                for (int i = 0; i < M; i++) {
                    Ioutj += (Vout(i*N + j) - Vout(i*N + j + M*N)) * G(i,j);
                }
                Iout.push_back(Ioutj);
            }

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