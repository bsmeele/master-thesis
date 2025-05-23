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

    crossbar.SetCrossbarParameters(3, INFINITY, INFINITY, 5, 3, 2);

    for (int it = 0; it < runs; it++) {
        std::vector<std::vector<bool>> weights(N, std::vector<bool>(M));
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                weights[i][j] = rand() % 2;  // rand() % 2 gives either 0 (false) or 1 (true)
            }
        }
        crossbar.SetRRAM(weights);

        if (print) {
            std::cout << "Weights:" << std::endl;
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < M; ++j) {
                    std::cout << weights[i][j] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

        float Vdd = 1.;
        Eigen::VectorXf Vappwl1 = Eigen::VectorXf::Random(M);
        Vappwl1 = (Vappwl1.array() > 0.5).select(Eigen::VectorXf::Constant(M, Vdd), Eigen::VectorXf::Zero(M));

        Eigen::VectorXf Vappbl2 = Eigen::VectorXf::Random(M);
        Vappbl2 = (Vappbl2.array() > 0.5).select(Eigen::VectorXf::Constant(M, Vdd), Eigen::VectorXf::Zero(M));
        // Eigen::VectorXf Vappbl2 = Eigen::VectorXf::Zero(M);

        Eigen::VectorXf Vappwl2 = Eigen::VectorXf::Zero(M);
        Eigen::VectorXf Vappbl1 = Eigen::VectorXf::Zero(M);

        if (print) {
            std::cout << "Vappwl1:\n" << Vappwl1 << std::endl << std::endl;
            std::cout << "Vappbl2:\n" << Vappbl2 << std::endl << std::endl;
        }

        Eigen::VectorXf V = Eigen::VectorXf::Zero(2*M*N);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                V(i*N + j) = Vappwl1(i);
            }
        }
        
        auto start_time = std::chrono::high_resolution_clock::now();

        Eigen::VectorXf Vout = crossbar.NonlinearSolve(V, Vappwl1, Vappwl2, Vappbl1, Vappbl2);

        auto end_time = std::chrono::high_resolution_clock::now();

        std::vector<float> Iout = crossbar.CalculateIout(Vout);

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

        auto execution_time = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
        std::cout << "Execution time: " << execution_time << " (us)" << std::endl;
        total_time += execution_time;

        if (print) {
            // std::cout << "Resistances:" << std::endl;
            // for (int i = 0; i < M; i++) {
            //     for (int j = 0; j < N; j++) {
            //         float v = Vout(i*N + j) - Vout(i*N + j + M*N);
            //         std::cout << crossbar.RRAM[i][j].GetResistance(v) << " ";
            //     }
            //     std::cout << std::endl;
            // }
            // std::cout << std::endl;

            std::cout << "Iout:" << std::endl;
            for (int j = 0; j < N; j++) {
                std::cout << Iout[j] << std::endl;
            }
            std::cout << std::endl;
        }
    }
    
    if (runs > 1) {
        std::cout << "Average execution time: " << total_time/runs << " us" << std::endl;
        // std::cout << "Average norm: " << total_norm/runs << std::endl;
        // std::cout << "Average iterations: " << (float) total_it/runs << std::endl;
    }



    // std::ofstream outfile("out.out");

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

    // // outfile << "t V I Nreal Treal Vschottky Vdiscplugserial Rschottky Rdisc Rplug Rseries Rtotal" << std::endl;

    // for (int i = 0; i < M; i++) {
    //     for (int j = 0; j < N; j++) {
    //         if (i == 0) {
    //             crossbar.access_transistors[i][j] = true;
    //         } else {
    //             crossbar.access_transistors[i][j] = false;
    //         }
    //         crossbar.RRAM[i][j].SetWeight(false);
    //     }
    // }

    // for (int i = 1; i < Vwave.size(); i++) {
    //     double dv = (Vwave[i][0] - Vappwl1(0)) / ((Vwave[i][1] - t) / dt);
    //     while (t < Vwave[i][1]) {
    //         // V = crossbar.NonlinearSolve(V, Vappwl1, Vappwl2, Vappbl1, Vappbl2);
    //         std::vector<std::vector<float>> I = crossbar.ApplyVoltage(V, Vappwl1, Vappwl2, Vappbl1, Vappbl2, dt);

    //         outfile << t;
    //         for (int n = 0; n < N; n++) {
    //             outfile << " " << crossbar.RRAM[0][n].Nreal;
    //         }
    //         outfile << std::endl;

    //         // float v = V(N-1) - V(N-1 + M*N);
    //         // outfile << t << " " << v << " " << I[0][N-1] << " " << crossbar.RRAM[0][N-1].Nreal << " " << crossbar.RRAM[0][N-1].Treal
    //         //     << " " << (v - (crossbar.RRAM[0][N-1].Rdisc + crossbar.RRAM[0][N-1].Rplug + crossbar.RRAM[0][N-1].Rseries) * I[0][N-1]) << " " << (crossbar.RRAM[0][N-1].Rdisc + crossbar.RRAM[0][N-1].Rplug + crossbar.RRAM[0][N-1].Rseries) * I[0][N-1]
    //         //     << " " << (v - (crossbar.RRAM[0][N-1].Rdisc + crossbar.RRAM[0][N-1].Rplug + crossbar.RRAM[0][N-1].Rseries) * I[0][N-1])/I[0][N-1] << " " << crossbar.RRAM[0][N-1].Rdisc << " " << crossbar.RRAM[0][N-1].Rplug << " " << crossbar.RRAM[0][N-1].Rseries
    //         //     << " " << v/I[0][N-1]
    //         //     << std::endl;

    //         // float v = V((M-1)*N) - V((M-1)*N + M*N);
    //         // outfile << t << " " << v << " " << I[M-1][0] << " " << crossbar.RRAM[M-1][0].Nreal << " " << crossbar.RRAM[M-1][0].Treal
    //         //     << " " << (v - (crossbar.RRAM[M-1][0].Rdisc + crossbar.RRAM[M-1][0].Rplug + crossbar.RRAM[M-1][0].Rseries) * I[M-1][0]) << " " << (crossbar.RRAM[M-1][0].Rdisc + crossbar.RRAM[M-1][0].Rplug + crossbar.RRAM[M-1][0].Rseries) * I[M-1][0]
    //         //     << " " << (v - (crossbar.RRAM[M-1][0].Rdisc + crossbar.RRAM[M-1][0].Rplug + crossbar.RRAM[M-1][0].Rseries) * I[M-1][0])/I[M-1][0] << " " << crossbar.RRAM[M-1][0].Rdisc << " " << crossbar.RRAM[M-1][0].Rplug << " " << crossbar.RRAM[M-1][0].Rseries
    //         //     << " " << v/I[M-1][0]
    //         //     << std::endl;

    //         Vappwl1(0) += dv;
    //         // for (int i = 0; i < N-1; i++) {
    //         //     Vappbl2(i) += dv;
    //         // }
    //         t += dt;
    //     }
    // }

    // outfile.close();
}
