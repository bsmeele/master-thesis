#include "nonlinear_crossbar_solver.h"
#include "crossbar_model/linear_crossbar_solver.h"
#include "crossbar_simulator.h"

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
    int total_it = 0;
    float total_norm = 0;

    CrossbarSimulator crossbar = CrossbarSimulator(M, N);

    for (int i = 0; i < runs; i++) {
        float Vdd = 1.;
        Eigen::VectorXf Vappwl1 = Eigen::VectorXf::Random(M);
        Vappwl1 = (Vappwl1.array() > 0.5).select(Eigen::VectorXf::Constant(M, Vdd), Eigen::VectorXf::Zero(M));

        Eigen::VectorXf Vappbl2 = Eigen::VectorXf::Random(M);
        Vappwl1 = (Vappbl2.array() > 0.5).select(Eigen::VectorXf::Constant(M, Vdd), Eigen::VectorXf::Zero(M));

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

        std::vector<float> Iout = crossbar.NonlinearSolve(V, Vappwl1, Vappwl2, Vappbl1, Vappbl2);

        auto end_time = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < Iout.size(); i++) {
                if (std::isnan(Iout[i])) {
                    std::cout << "NAN detected" << std::endl;
                    assert(false);
                }
                if (std::isinf(Iout[i])) {
                    std::cout << "Inf detected" << std::endl;
                    assert(false);
                }
        }

        auto execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        std::cout << "Execution time: " << execution_time << " (ms)" << std::endl;
        total_time += execution_time;

        if (print) {
            std::cout << "Iout:" << std::endl;
            for (int j = 0; j < N; j++) {
                std::cout << Iout[j] << std::endl;
            }
            std::cout << std::endl;
        }
    }
    
    if (runs > 1) {
        std::cout << "Average execution time: " << total_time/runs << " ms" << std::endl;
        // std::cout << "Average norm: " << total_norm/runs << std::endl;
        // std::cout << "Average iterations: " << (float) total_it/runs << std::endl;
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

    // Vappwl1(0) = Vwave[0][0];
    // double t = Vwave[0][1];

    // Eigen::VectorXf V = Eigen::VectorXf::Zero(2*M*N);

    // outfile << "t V I Nreal Treal Vschottky Vdiscplugserial Rschottky Rdisc Rplug Rseries Rtotal" << std::endl;

    // for (int i = 0; i < M; i++) {
    //     for (int j = 0; j < N; j++) {
    //         access_transistors[i][j] = false;
    //     }
    // }

    // access_transistors[0][N-1] = true;

    // for (int i = 1; i < Vwave.size(); i++) {
    //     double dv = (Vwave[i][0] - Vappwl1(0)) / ((Vwave[i][1] - t) / dt);
    //     while (t < Vwave[i][1]) {
    //         V = FixedpointSolve(RRAM, access_transistors, V, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, print);

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

    //                 if (i == 0 && j == N-1) {
    //                     outfile << t << " " << v << " " << I << " " << RRAM[i][j].Nreal << " " << RRAM[i][j].Treal
    //                     << " " << (v - (RRAM[i][j].Rdisc + RRAM[i][j].Rplug + RRAM[i][j].Rseries) * I) << " " << (RRAM[i][j].Rdisc + RRAM[i][j].Rplug + RRAM[i][j].Rseries) * I
    //                     << " " << (v - (RRAM[i][j].Rdisc + RRAM[i][j].Rplug + RRAM[i][j].Rseries) * I)/I << " " << RRAM[i][j].Rdisc << " " << RRAM[i][j].Rplug << " " << RRAM[i][j].Rseries
    //                     << " " << v/I
    //                     << std::endl;
    //                 }
    //             }
    //         }

    //         Vappwl1(0) += dv;
    //         t += dt;
    //     }
    // }

    // outfile.close();
}
