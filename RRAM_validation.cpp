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

    int M = (argc >= 2) ? std::atoi(argv[1]) : 32;
    int N = (argc >= 3) ? std::atoi(argv[2]) : 32;

    bool print = (argc >= 4 && std::atoi(argv[3]) == 1) ? true : false;

    // 1 ohm resistors
    // Ndiscmin of 0.0001, default parameters otherwise
    // Transistor gate is connected to source line, thus memristor will be connected if source line is 1
    CrossbarSimulator crossbar = CrossbarSimulator(M, N);
    for (int i = 0; i < crossbar.RRAM.size(); i++) {
        for (int j = 0; j < crossbar.RRAM[i].size(); j++) {
            crossbar.RRAM[i][j].Ndiscmin = 0.0001;
            crossbar.RRAM[i][j].Ninit = 0.0001;
        }
    }

    std::vector<std::vector<bool>> weights;
    for (int i = 0; i < M; i++) {
        std::vector<bool> row;
        for (int j = 0; j < N; j++) {
            if (std::rand()%2 == 0) {
                row.push_back(true);
            } else {
                row.push_back(false);
            }
        }
        weights.push_back(row);
    }
    crossbar.SetRRAM(weights);

    std::vector<bool> Vapp;
    for (int i = 0; i < M; i++) {
        if (std::rand()%2 == 0) {
            Vapp.push_back(true);
        } else {
            Vapp.push_back(false);
            crossbar.access_transistors[i] = std::vector<bool>(N, false);
        }
    }
    Vapp[0] = true;

    // 0.1 V read voltage
    // width of 50 us
    // 5 us rise and fall time
    float Vdd = 0.1;
    Eigen::VectorXf Vappwl1 = Eigen::VectorXf::Zero(M);
    Eigen::VectorXf Vappwl2 = Eigen::VectorXf::Zero(M);
    Eigen::VectorXf Vappbl1 = Eigen::VectorXf::Zero(M);
    Eigen::VectorXf Vappbl2 = Eigen::VectorXf::Zero(M);

    Eigen::VectorXf Vguess = Eigen::VectorXf::Zero(2*M*N);

    float dt = 1e-6;
    std::vector<std::array<float, 2>> Vwave = { {{0, 0}}, {{Vdd, 5e-6}}, {{Vdd, 45e-6}}, {{0, 50e-6}}};

    float V = 0;
    float t = 0;

    std::vector<std::vector<std::vector<float>>> Iwave;

    auto start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < Vwave.size(); i++) {
        float dv = (Vwave[i][0] - V) / ((Vwave[i][1] - t) / dt);
        while (t < Vwave[i][1]) {
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < N; j++) {
                    Vguess(i*N + j) = Vappwl1(i);
                }
            }

            std::vector<std::vector<float>> Iout = crossbar.ApplyVoltage(Vguess, Vappwl1, Vappwl2, Vappbl1, Vappbl2, dt);
            Iwave.push_back(Iout);

            for (int j = 0; j < M; j++) {
                if (Vapp[j]) {
                    Vappwl1(j) += dv;
                }
            }
            V += dv;
            t += dt;
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "Execution time: " << execution_time << " (ms)" << std::endl;

    std::vector<float> Iout;
    for (int n = 0; n < N; n++) {
        float Iavg = 0;
        for (int i = 0; i < 40; i++) {
            for (int m = 0; m < M; m++) {
                Iavg += Iwave[i][m][n];
            }
        }
        Iavg /= 40.;
        Iout.push_back(Iavg);
    }

    if (print) {
        std::cout << "Vapplied (boolean):" << std::endl;
        for (int i = 0; i < Vapp.size(); i++) {
            std::cout << Vapp[i] << std::endl;
        }
        std::cout << std::endl;

        std::cout << "Weights (boolean):" << std::endl;
        for (int i = 0; i < weights.size(); i++) {
            for (int j = 0; j < weights[0].size(); j++) {
                std::cout << weights[i][j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;

        std::cout << "Iout:" << std::endl;
        for (int i = 0; i < Iout.size(); i++) {
            std::cout << Iout[i] << std::endl;
        }
    }
}
