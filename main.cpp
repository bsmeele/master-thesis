#include "RRAM_crossbar_model.h"
#include "crossbar_model/crossbar_array_model.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <iostream>
#include <chrono>

int main(int argc, char* argv[]) {
    srand((unsigned int) time(0));

    int M = (argc >= 2) ? std::atoi(argv[1]) : 16;
    int N = (argc >= 3) ? std::atoi(argv[2]) : 16;

    int runs = (argc >= 4) ? std::atoi(argv[3]) : 1;

    bool print = (argc >= 5 && std::atoi(argv[4]) == 1) ? true : false;

    std::cout << "Solving " << M << "x" << N << " system in " << runs << " run(s)" << std::endl;

    long long total_time = 0;

    float Rswl1 = 3.;
    float Rswl2 = INFINITY;
    float Rsbl1 = INFINITY;
    float Rsbl2 = 5.;

    float Rwl = 3.;
    float Rbl = 2.;

    Eigen::SparseMatrix<float> G_ABCD = PartiallyPrecomputeG_ABCD(M, N, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);

    std::vector<std::vector<JART_VCM_v1b_var>> RRAM;
    for (int i = 0; i < M; i++) {
        std::vector<JART_VCM_v1b_var> row;
        for (int j = 0; j < N; j++) {
            row.push_back(JART_VCM_v1b_var());
        }
        RRAM.push_back(row);
    }

    for (int i = 0; i < runs; i++) {
        float Vdd = 1.;
        Eigen::VectorXf Vappwl1 = Eigen::VectorXf::Random(M);
        Vappwl1 = (Vappwl1.array() > 0.5).select(Eigen::VectorXf::Constant(M, Vdd), Eigen::VectorXf::Zero(M));

        if (runs == 1 && print && M == 16 && N == 16) {
            Vappwl1 = Eigen::VectorXf::Zero(M);
            Vappwl1(0) = 1.;
            Vappwl1(12) = 1.;
            Vappwl1(13) = 1.;
        }

        Eigen::VectorXf Vappwl2 = Eigen::VectorXf::Zero(M);
        Eigen::VectorXf Vappbl1 = Eigen::VectorXf::Zero(M);
        Eigen::VectorXf Vappbl2 = Eigen::VectorXf::Zero(M);

        if (print) {
            std::cout << "Vappwl1:\n" << Vappwl1 << std::endl << std::endl;
        }

        Eigen::VectorXf V = Eigen::VectorXf::Zero(2*M*N);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                V(i*N + j) = Vappwl1(i);
            }
        }
        auto start_time = std::chrono::high_resolution_clock::now();

        Eigen::VectorXf Vout = FixedpointSolve(RRAM, V, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, print);
        // Eigen::VectorXf Vout = NewtonRaphsonSolve(RRAM, V, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, print);
        // Eigen::VectorXf Vout = BroydenSolve(RRAM, V, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, print);
        // Eigen::VectorXf Vout = BroydenInvSolve(RRAM, V, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, print);

        auto end_time = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < Vout.size(); i++) {
            if (std::isnan(Vout(i))) {
                std::cout << "NAN detected" << std::endl;
                assert(false);
            }
            if (std::isinf(Vout(i))) {
                std::cout << "Inf detected" << std::endl;
                assert(false);
            }
        }

        auto execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        std::cout << "Execution time: " << execution_time << " (ms)" << std::endl;
        total_time += execution_time;

        if (print) {
            Eigen::MatrixXf G(M, N);
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < N; j++) {
                    float v = Vout(i*N + j) - Vout(i*N + j + M*N);
                    G(i, j) = 1./RRAM[i][j].GetResistance(v);
                }
            }

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
    }
    
    if (runs > 1) {
        std::cout << "Average execution time: " << total_time/runs << " ms" << std::endl;
        // std::cout << "Average norm: " << avg_norm/runs << std::endl;
        // std::cout << "Average iterations: " << (float) avg_it/runs << std::endl;
    }
}