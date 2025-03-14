#include "linear_crossbar_solver.h"

#include <iostream>
#include <chrono>
#include <cassert>
#include <string>
#include <sstream>
#include <fstream>

void LoadFromFile(
    std::string filepath,
    int& M, int& N,
    Eigen::MatrixXf& G,
    Eigen::VectorXf& Vguess, Eigen::SparseMatrix<float>& G_ABCD,
    Eigen::VectorXf& Vappwl1, Eigen::VectorXf& Vappwl2,
    Eigen::VectorXf& Vappbl1, Eigen::VectorXf& Vappbl2,
    float& Rswl1, float& Rswl2, float& Rsbl1, float& Rsbl2,
    float& Rwl, float& Rbl
) {
    std::ifstream file(filepath);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open the file." << std::endl;
        assert(false);
    }

    G.setZero();
    Vguess.setZero();
    G_ABCD.setZero();
    Vappwl1.setZero();
    Vappwl2.setZero();
    Vappbl1.setZero();
    Vappbl2.setZero();


    int flag = 0;
    int flag2 = 0;
    std::string line;
    std::stringstream stream;
    while (std::getline(file, line)) {
        if (line.empty()) { continue; }
        switch (flag) {
            case 0:
                M = std::stoi(line);
                flag = 1;
                break;
            case 1:
                N = std::stoi(line);
                flag = 2;
                break;
            case 2:
                stream.clear();
                stream.str(line);
                for (int i = 0; i < N; i++) {
                    float g;
                    stream >> g;
                    G(flag2, i) = g;
                }
                flag2 += 1;
                if (flag2 == M) {
                    flag = 3;
                    flag2 = 0;
                }
                break;
            case 3:
                Vguess(flag2) = std::stof(line);
                flag2 += 1;
                if (flag2 == 2*M*N) {
                    flag = 4;
                    flag2 = 0;
                }
                break;
            case 4:
                stream.clear();
                stream.str(line);
                for (int i = 0; i < 2*M*N; i++) {
                    float value;
                    stream >> value;
                    if (value != 0) {
                        G_ABCD.insert(flag2, i) = value;
                    }
                }
                flag2 += 1;
                if (flag2 == 2*M*N) {
                    flag = 5;
                    flag2 = 0;
                }
                break;
            case 5:
                Vappwl1(flag2) = std::stof(line);
                flag2 += 1;
                if (flag2 == M) {
                    flag = 6;
                    flag2 = 0;
                }
                break;
            case 6:
                Vappwl2(flag2) = std::stof(line);
                flag2 += 1;
                if (flag2 == M) {
                    flag = 7;
                    flag2 = 0;
                }
                break;
            case 7:
                Vappbl1(flag2) = std::stof(line);
                flag2 += 1;
                if (flag2 == N) {
                    flag = 8;
                    flag2 = 0;
                }
                break;
            case 8:
                Vappbl2(flag2) = std::stof(line);
                flag2 += 1;
                if (flag2 == N) {
                    flag = 9;
                    flag2 = 0;
                }
                break;
            case 9:
                Rswl1 = std::stof(line);
                flag = 10;
                break;
            case 10:
                Rswl2 = std::stof(line);
                flag = 11;
                break;
            case 11:
                Rsbl1 = std::stof(line);
                flag = 12;
                break;
            case 12:
                Rsbl2 = std::stof(line);
                flag = 13;
                break;
            case 13:
                Rwl = std::stof(line);
                flag = 14;
                break;
            case 14:
                Rbl = std::stof(line);
                flag = 15;
                break;
        }
    }
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

    Eigen::VectorXf Vref = Eigen::VectorXf::Zero(2*M*N);

    float Rswl1 = 3.;
    float Rswl2 = INFINITY;
    float Rsbl1 = INFINITY;
    float Rsbl2 = 5.;

    float Rwl = 3.;
    float Rbl = 2.;

    Eigen::SparseMatrix<float> G_ABCD = PartiallyPrecomputeG_ABCD(M, N, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);

    for (int i = 0; i < runs; i++) {
        float Rmin = 1000.;
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
        // if (M == 3 && N == 3) {
        //     R = Eigen::MatrixXf::Zero(M, N);
        //     // R(0, 0) = 10.;
        //     // R(0, 1) = 15.;
        //     // R(0, 2) = 20.;
        //     // R(1, 0) = 25.;
        //     // R(1, 1) = 30.;
        //     // R(1, 2) = 35.;
        //     // R(2, 0) = 40.;
        //     // R(2, 1) = 45.;
        //     // R(2, 2) = 50.;
    
        //     R(0, 0) = 1e4;
        //     R(0, 1) = 2e4;
        //     R(0, 2) = 3e4;
        //     R(1, 0) = 4e4;
        //     R(1, 1) = 5e4;
        //     R(1, 2) = 6e4;
        //     R(2, 0) = 7e4;
        //     R(2, 1) = 8e4;
        //     R(2, 2) = 9e4;
        // }
        // R(0, 0) = 1e4;
        // R(0, 1) = 2e4;
        // R(0, 2) = 3e4;
        // R(0, 3) = 4e4;
        // R(0, 4) = 5e4;
        // R(0, 5) = 6e4;
        // R(0, 6) = 7e4;
        // R(0, 7) = 8e4;
        // R(1, 0) = 9e4;
        // R(1, 1) = 10e4;
        // R(1, 2) = 11e4;
        // R(1, 3) = 12e4;
        // R(1, 4) = 13e4;
        // R(1, 5) = 14e4;
        // R(1, 6) = 15e4;
        // R(1, 7) = 16e4;
        // R(2, 0) = 17e4;
        // R(2, 1) = 18e4;
        // R(2, 2) = 19e4;
        // R(2, 3) = 20e4;
        // R(2, 4) = 21e4;
        // R(2, 5) = 22e4;
        // R(2, 6) = 23e4;
        // R(2, 7) = 24e4;
        // R(3, 0) = 25e4;
        // R(3, 1) = 26e4;
        // R(3, 2) = 27e4;
        // R(3, 3) = 28e4;
        // R(3, 4) = 29e4;
        // R(3, 5) = 30e4;
        // R(3, 6) = 31e4;
        // R(3, 7) = 32e4;
        // R(4, 0) = 33e4;
        // R(4, 1) = 34e4;
        // R(4, 2) = 35e4;
        // R(4, 3) = 36e4;
        // R(4, 4) = 37e4;
        // R(4, 5) = 38e4;
        // R(4, 6) = 39e4;
        // R(4, 7) = 40e4;
        // R(5, 0) = 41e4;
        // R(5, 1) = 42e4;
        // R(5, 2) = 43e4;
        // R(5, 3) = 44e4;
        // R(5, 4) = 45e4;
        // R(5, 5) = 46e4;
        // R(5, 6) = 47e4;
        // R(5, 7) = 48e4;
        // R(6, 0) = 49e4;
        // R(6, 1) = 50e4;
        // R(6, 2) = 51e4;
        // R(6, 3) = 52e4;
        // R(6, 4) = 53e4;
        // R(6, 5) = 54e4;
        // R(6, 6) = 55e4;
        // R(6, 7) = 56e4;
        // R(7, 0) = 57e4;
        // R(7, 1) = 58e4;
        // R(7, 2) = 59e4;
        // R(7, 3) = 60e4;
        // R(7, 4) = 61e4;
        // R(7, 5) = 62e4;
        // R(7, 6) = 63e4;
        // R(7, 7) = 64e4;

        if (print) {
            std::cout << "R:\n" << R << std::endl << std::endl;
        }

        Eigen::MatrixXf G = R.cwiseInverse();

        if (print) {
            std::cout << "G:\n" << G << std::endl << std::endl;
        }

        Eigen::VectorXf Vappwl1 = Eigen::VectorXf::Random(M);
        Vappwl1 = (Vappwl1.array() > 0.5).select(Eigen::VectorXf::Constant(M, Vdd), Eigen::VectorXf::Zero(M));
        // if (M == 3 && N == 3) {
        //     Vappwl1(0) = 0.5;
        //     Vappwl1(1) = 1;
        //     Vappwl1(2) = 1.5;
        // }
        // Vappwl1(0) = 0.5;
        // Vappwl1(1) = 1;
        // Vappwl1(2) = 1.5;
        // Vappwl1(3) = 1;
        // Vappwl1(4) = 0.5;
        // Vappwl1(5) = 1;
        // Vappwl1(6) = 1.5;
        // Vappwl1(7) = 1;

        Eigen::VectorXf Vappwl2 = Eigen::VectorXf::Zero(M);
        Eigen::VectorXf Vappbl1 = Eigen::VectorXf::Zero(M);
        Eigen::VectorXf Vappbl2 = Eigen::VectorXf::Zero(M);

        if (print) {
            std::cout << "Vappwl1:\n" << Vappwl1 << std::endl << std::endl;
        }

        // Set the initial guess to be the wordline voltage
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                V(i*N + j) = Vappwl1(i);
            }
        }

        // if (M == 16 && N == 16 && print) {
        //     load_from_file("../out.txt", M, N, G, V, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);
        //     std::cout << "Inputs loaded from file" << std::endl << std::endl;
        // }

        auto start_time = std::chrono::high_resolution_clock::now();

        Eigen::VectorXf Vout = SolveCam(G, V, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, false);

        auto end_time = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < Vout.size(); i++) {
            assert(!std::isnan(Vout(i)));
        }

        if (print) {
            std::cout << "Vout:\n" << Vout << std::endl << std::endl;

            std::cout << "Norm: " << (V - Vout).norm() << std::endl << std::endl;
            
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