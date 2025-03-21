#include "nonlinear_crossbar_solver.h"
#include "crossbar_model/linear_crossbar_solver.h"
#include "crossbar_simulator.h"
#include "simulation_settings.h"

#include </home/earapidis/dev/eigen/Eigen/Dense>
#include </home/earapidis/dev/eigen/Eigen/Sparse>

#include <iostream>
#include <chrono>
#include <fstream>
#include <filesystem>
#include <cstdint>
#include <algorithm>

namespace fs = std::filesystem;

std::vector<std::vector<int64_t>> readMatrixFromFile(const std::filesystem::path& filePath) {
    std::ifstream file(filePath, std::ios::binary);
    
    if (!file) {
        throw std::runtime_error("Failed to open file: " + filePath.string());
    }

    // Read rows and columns
    int64_t rows, cols;
    file.read(reinterpret_cast<char*>(&rows), sizeof(int64_t));
    file.read(reinterpret_cast<char*>(&cols), sizeof(int64_t));

    // Read data
    std::vector<std::vector<int64_t>> matrix(rows, std::vector<int64_t>(cols));

    // Read the data into the matrix
    for (int64_t i = 0; i < rows; ++i) {
        file.read(reinterpret_cast<char*>(matrix[i].data()), cols * sizeof(int64_t));
    }

    return matrix;
}



int main(int argc, char* argv[]) {
    // srand((unsigned int) time(0));

    // int M = (argc >= 2) ? std::atoi(argv[1]) : 32;
    // int N = (argc >= 3) ? std::atoi(argv[2]) : 32;

    // bool print = (argc >= 4 && std::atoi(argv[3]) == 1) ? true : false;

    std::filesystem::path top_dir, input_path, weight_path, mac_path, mem_path;
    if (argc >= 2) {
        top_dir = fs::absolute(argv[1]);
    }

    input_path = top_dir/"input.bin";
    weight_path = top_dir/"weight.bin";
    mac_path = top_dir/"mac.bin";
    mem_path = top_dir/"mem.bin";


    auto input_data = readMatrixFromFile(input_path);
    auto weight_data = readMatrixFromFile(weight_path);

    int M = weight_data.size();
    int N = weight_data[0].size();

    CrossbarSimulator crossbar = CrossbarSimulator(M, N);

    std::vector<std::vector<bool>> weights;
    for (int i = 0; i < M; i++) {
        std::vector<bool> row;
        for (int j = 0; j < N; j++) {
            if (weight_data[i][j] != 0) {
                row.push_back(true);
            } else {
                row.push_back(false);
            }
        }
        weights.push_back(row);
    }

    std::vector<std::vector<std::vector<float>>> output_data;
    std::vector<std::vector<float>> output_data_MAC;

    auto start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < input_data.size(); i++) {
        crossbar.SetRRAM(weights);

        std::vector<bool> Vapp;
        for (int j = 0; j < M; j++) {
            if (input_data[i][j] != 0) {
                Vapp.push_back(true);
                crossbar.access_transistors[j] = std::vector<bool>(N, true);
            } else {
                Vapp.push_back(false);
                crossbar.access_transistors[j] = std::vector<bool>(N, false);
            }
        }

        float Vdd = 0.1;
        Eigen::VectorXf Vappwl1 = Eigen::VectorXf::Zero(M);
        Eigen::VectorXf Vappwl2 = Eigen::VectorXf::Zero(M);
        Eigen::VectorXf Vappbl1 = Eigen::VectorXf::Zero(M);
        Eigen::VectorXf Vappbl2 = Eigen::VectorXf::Zero(M);

        Eigen::VectorXf Vguess = Eigen::VectorXf::Zero(2*M*N);

        float dt = simulation_time_step;
        std::vector<std::array<float, 2>> Vwave = {
            {{0, 0}},
            {{voltage_pulse_height, voltage_pulse_rise_time}},
            {{voltage_pulse_height, voltage_pulse_width - voltage_pulse_fall_time}},
            {{0, voltage_pulse_width}}
        };

        float V = 0;
        float t = 0;

        std::vector<std::vector<std::vector<float>>> Iwave;

        for (int i = 0; i < Vwave.size(); i++) {
            float dv = (Vwave[i][0] - V) / ((Vwave[i][1] - t) / dt);
            while (t < Vwave[i][1]) {
                for (int i = 0; i < M; i++) {
                    for (int j = 0; j < N; j++) {
                        Vguess(i*N + j) = Vappwl1(i);
                    }
                }

                Iwave.push_back(crossbar.ApplyVoltage(Vguess, Vappwl1, Vappwl2, Vappbl1, Vappbl2, dt));

                for (int j = 0; j < M; j++) {
                    if (Vapp[j]) {
                        Vappwl1(j) += dv;
                    }
                }
                V += dv;
                t += dt;
            }
        }

        std::vector<std::vector<float>> Iout_avg;
        for (int m = 0; m < M; m++) {
            std::vector<float> row;
            for (int n = 0; n < N; n++) {
                float Iavg = 0;
                for (int j = voltage_pulse_rise_time/simulation_time_step; j < (voltage_pulse_width - voltage_pulse_fall_time)/simulation_time_step; j++) {
                    Iavg += Iwave[j][m][n];
                }
                Iavg /= (voltage_pulse_width - voltage_pulse_rise_time - voltage_pulse_fall_time) / simulation_time_step;
                row.push_back(Iavg);
            }
            Iout_avg.push_back(row);
        }
        // std::vector<std::vector<float>> Iout_avg = Iwave[25];

        std::vector<float> Iout_MAC;
        for (int n = 0; n < N; n++) {
            float IMAC = 0;
            for (int m = 0; m < M; m++) {
                IMAC += Iout_avg[m][n];
            }
            Iout_MAC.push_back(IMAC);
        }

        output_data.push_back(Iout_avg);
        output_data_MAC.push_back(Iout_MAC);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "Execution time: " << execution_time << " (ms)" << std::endl;


    std::ofstream outfile(top_dir/"output.bin", std::ios::binary);

    if (!outfile) {
        std::cout << "Failed to open file: " << top_dir/"output.bin" << std::endl;
        return 1;
    }

    int64_t depth = output_data.size();
    int64_t rows = output_data.empty() ? 0 : output_data[0].size();
    int64_t cols = rows ? output_data[0][0].size() : 0;

    outfile.write(reinterpret_cast<const char*>(&depth), sizeof(int64_t));
    outfile.write(reinterpret_cast<const char*>(&rows), sizeof(int64_t));
    outfile.write(reinterpret_cast<const char*>(&cols), sizeof(int64_t));

    for (int64_t i = 0; i < depth; ++i) {
        for (int64_t j = 0; j < rows; ++j) {
            outfile.write(reinterpret_cast<const char*>(output_data[i][j].data()), cols * sizeof(float));
        }
    }

    outfile.close();

    outfile = std::ofstream(top_dir/"out_mac.bin", std::ios::binary);

    if (!outfile) {
        std::cout << "Failed to open file: " << top_dir/"MAC.bin" << std::endl;
        return 1;
    }

    rows = output_data_MAC.size();
    cols = output_data_MAC.empty() ? 0 : output_data_MAC[0].size();

    outfile.write(reinterpret_cast<const char*>(&rows), sizeof(int64_t));
    outfile.write(reinterpret_cast<const char*>(&cols), sizeof(int64_t));

    for (int64_t i = 0; i < rows; ++i) {
        outfile.write(reinterpret_cast<const char*>(output_data_MAC[i].data()), cols * sizeof(float));
    }


    outfile.close();
    
}
