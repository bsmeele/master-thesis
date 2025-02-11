#include "RRAM_crossbar_model.h"
#include "crossbar_model/crossbar_array_model.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <iostream>
#include <chrono>
#include <fstream>

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

        Eigen::VectorXf Vappbl2 = Eigen::VectorXf::Random(M);
        Vappwl1 = (Vappbl2.array() > 0.5).select(Eigen::VectorXf::Constant(M, Vdd), Eigen::VectorXf::Zero(M));

        // if (runs == 1 && print && M == 16 && N == 16) {
        //     Vappwl1 = Eigen::VectorXf::Zero(M);
        //     Vappbl2 = Eigen::VectorXf::Zero(M);
        // }

        Eigen::VectorXf Vappwl2 = Eigen::VectorXf::Zero(M);
        Eigen::VectorXf Vappbl1 = Eigen::VectorXf::Zero(M);

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



    // std::ofstream outfile("out.txt");

    // if (!outfile) {
    //     std::cout << "No out file" << std::endl;
    //     return 1;
    // }

    // const double dt = 1e-3;
    
    // Eigen::VectorXf Vappwl1 = Eigen::VectorXf::Zero(M);
    // Eigen::VectorXf Vappwl2 = Eigen::VectorXf::Zero(M);
    // Eigen::VectorXf Vappbl1 = Eigen::VectorXf::Zero(M);
    // Eigen::VectorXf Vappbl2 = Eigen::VectorXf::Zero(M);

    // std::vector<std::array<double, 2>> Vwave;
    // Vwave.push_back({0, 0});
    // Vwave.push_back({-1.5, 1.5});
    // Vwave.push_back({0, 3});
    // Vwave.push_back({1.5, 4.5});
    // Vwave.push_back({0, 6});

    // Vappwl1(M-1) = Vwave[0][0];
    // double t = Vwave[0][1];

    // Eigen::VectorXf V = Eigen::VectorXf::Zero(2*M*N);

    // outfile << "t V I Nreal Treal Vschottky Vdiscplugserial Rschottky Rdisc Rplug Rseries Rtotal" << std::endl;

    // for (int i = 1; i < Vwave.size(); i++) {
    //     double dv = (Vwave[i][0] - Vappwl1(M-1)) / ((Vwave[i][1] - t) / dt);
    //     while (t < Vwave[i][1]) {
    //         V = FixedpointSolve(RRAM, V, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, print);

    //         for (int i = 0; i < V.size(); i++) {
    //             if (std::isnan(V(i))) {
    //                 std::cout << "NAN detected" << std::endl;
    //                 assert(false);
    //             }
    //             if (std::isinf(V(i))) {
    //                 std::cout << "Inf detected" << std::endl;
    //                 assert(false);
    //             }
    //         }

    //         for (int i = 0; i < M; i++) {
    //             for (int j = 0; j < N; j++) {
    //                 float v = V(i*N + j) - V(i*N + j + M*N);
    //                 double I = RRAM[i][j].ApplyVoltage(v, dt);

    //                 if (std::isnan(I)) {
    //                     std::cout << "NAN detected" << std::endl;
    //                     assert(false);
    //                 }
    //                 if (std::isinf(I)) {
    //                     std::cout << "Inf detected" << std::endl;
    //                     assert(false);
    //                 }

    //                 if (i == M-1 && j == 0) {
    //                     outfile << t << " " << v << " " << I << " " << RRAM[i][j].Nreal << " " << RRAM[i][j].Treal
    //                     << " " << (v - (RRAM[i][j].Rdisc + RRAM[i][j].Rplug + RRAM[i][j].Rseries) * I) << " " << (RRAM[i][j].Rdisc + RRAM[i][j].Rplug + RRAM[i][j].Rseries) * I
    //                     << " " << (v - (RRAM[i][j].Rdisc + RRAM[i][j].Rplug + RRAM[i][j].Rseries) * I)/I << " " << RRAM[i][j].Rdisc << " " << RRAM[i][j].Rplug << " " << RRAM[i][j].Rseries
    //                     << " " << v/I
    //                     << std::endl;
    //                 }
    //             }
    //         }

    //         Vappwl1(M-1) += dv;
    //         t += dt;
    //     }
    // }

    // outfile.close();
}