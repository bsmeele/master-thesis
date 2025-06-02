#include "nonlinear_crossbar_solver.h"
#include "crossbar_model/linear_crossbar_solver.h"
#include "crossbar_simulator.h"
#include "simulation_settings.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <iostream>
#include <chrono>
#include <fstream>
#include <filesystem>
#include <cstdint>
#include <algorithm>
#include <regex>
#include <string>

// Input argumets: M N batch_size write weight_ratio input_ratio
int main(int argc, char* argv[]) {
    srand((unsigned int) time(0));

    int M = (argc >= 2) ? std::atoi(argv[1]) : 64;
    int N = (argc >= 3) ? std::atoi(argv[2]) : 64;

    int batch_size = (argc >= 4) ? std::atoi(argv[3]) : 1;

    bool write = (argc >= 5 && std::atoi(argv[4]) == 1) ? true : false;

    float weight_ratio = (argc >= 6) ? std::atof(argv[5]) : 0.5;  // 0 is all HRS, 1 is all LRS
    float input_ratio = (argc >= 7) ? std::atof(argv[6]) : 0.5;  // 0 is all off, 1 is all on

    std::cout << "Simulating " << batch_size << " " << M << "X" << N << " crossbars" << std::endl;

    std::filesystem::path top_dir = "RRAM_test_data";
    if (write) {
        // File organization:
        // RRAM_test_data
        //   batch_n
        //     log.txt
        //     sim_n
        //       weight.bin
        //       input.bin
        //       out.bin
        //       out_MAC.bin

        std::regex batch_regex(R"(batch_(\d+))");
        int max_batch = -1;

        for (const auto& entry : std::filesystem::directory_iterator(top_dir)) {
            if (entry.is_directory()) {
                std::smatch match;
                std::string folder_name = entry.path().filename().string();

                if (std::regex_match(folder_name, match, batch_regex)) {
                    int n = std::stoi(match[1]);
                    if (n > max_batch) {
                        max_batch = n;
                    }
                }
            }
        }

        std::string batch_dir = "batch_" + std::to_string(max_batch + 1);
        top_dir = top_dir/batch_dir;
        std::filesystem::create_directory(top_dir);
    }
    
    // Initialize crossbar
    CrossbarSimulator crossbar = CrossbarSimulator(M, N);
    
    std::vector<std::vector<bool>> weights(M, std::vector<bool>(N));
    std::vector<bool> input(M);

    float total_time = 0.;

    float dt = simulation_time_step;
    std::vector<std::array<float, 2>> Vwave = {
        {{0, 0}},
        {{voltage_pulse_height, voltage_pulse_rise_time}},
        {{voltage_pulse_height, voltage_pulse_width - voltage_pulse_fall_time}},
        {{0, voltage_pulse_width}}
    };

    // for batch_num
    for (int i = 0; i < batch_size; i++) {
        std::ostringstream ss;
        ss << i << " out of " << batch_size << " sims completed (" << std::fixed << std::setprecision(1) << (float) i/batch_size * 100. << "%)";
        std::cout << "\r" << ss.str() << std::flush;
        
        // Randomly initialize the weigths based on the ratio
        for (int m = 0; m < M; m++) {
            for (int n = 0; n < N; n++) {
                weights[m][n] = ((double) rand())/RAND_MAX < weight_ratio;
            }
        }
        // Set weights
        crossbar.SetRRAM(weights);

        // Randomly initialize the input based on the ratio
        for (int m = 0; m < M; m++) {
            input[m] = ((double) rand())/RAND_MAX < input_ratio;
        }

        // Set access transistors based on the input voltage
        for (int m = 0; m < M; m++) {
            if (input[m]) {
                crossbar.access_transistors[m] = std::vector<bool>(N, true);
            } else {
                crossbar.access_transistors[m] = std::vector<bool>(N, false);
            }
        }
        
        Eigen::VectorXf Vappwl1 = Eigen::VectorXf::Zero(M);
        Eigen::VectorXf Vappwl2 = Eigen::VectorXf::Zero(M);
        Eigen::VectorXf Vappbl1 = Eigen::VectorXf::Zero(M);
        Eigen::VectorXf Vappbl2 = Eigen::VectorXf::Zero(M);

        Eigen::VectorXf Vguess = Eigen::VectorXf::Zero(2*M*N);

        float V = 0;
        float t = 0;

        std::vector<std::vector<std::vector<float>>> Iwave;

        auto start_time = std::chrono::high_resolution_clock::now();

        // Simulate
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
                    if (input[j]) {
                        Vappwl1(j) += dv;
                    }
                }

                V += dv;
                t += dt;
            }
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        total_time += execution_time;

        // Collect results
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

        std::vector<float> Iout_MAC;
        for (int n = 0; n < N; n++) {
            float IMAC = 0;
            for (int m = 0; m < M; m++) {
                IMAC += Iout_avg[m][n];
            }
            Iout_MAC.push_back(IMAC);
        }

        // Write all to files
        if (write) {
            std::string sim_dir = "sim_" + std::to_string(i);
            std::filesystem::create_directory(top_dir/sim_dir);

            std::ofstream weight_file(top_dir/sim_dir/"weight.bin");
            weight_file.write(reinterpret_cast<const char*>(&M), sizeof(int64_t));
            weight_file.write(reinterpret_cast<const char*>(&N), sizeof(int64_t));
            for (int m = 0; m < M; m++) {
                for (int n = 0; n < N; n++) {
                    int64_t value = weights[m][n] ? 1 : 0;
                    weight_file.write(reinterpret_cast<const char*>(&value), sizeof(int64_t));
                }
            }
            weight_file.close();

            std::ofstream input_file(top_dir/sim_dir/"input.bin");
            input_file.write(reinterpret_cast<const char*>(&M), sizeof(int64_t));
            for (int m = 0; m < M; m++) {
                int64_t value = input[m] ? 1 : 0;
                input_file.write(reinterpret_cast<const char*>(&value), sizeof(int64_t));
            }
            input_file.close();

            std::ofstream out_file(top_dir/sim_dir/"out.bin");
            out_file.write(reinterpret_cast<const char*>(&M), sizeof(int64_t));
            out_file.write(reinterpret_cast<const char*>(&N), sizeof(int64_t));
            for (int m = 0; m < M; m++) {
                out_file.write(reinterpret_cast<const char*>(Iout_avg[m].data()), N * sizeof(float));
            }
            out_file.close();
            
            std::ofstream mac_file(top_dir/sim_dir/"out_MAC.bin");
            mac_file.write(reinterpret_cast<const char*>(&M), sizeof(int64_t));
            mac_file.write(reinterpret_cast<const char*>(Iout_MAC.data()), N * sizeof(float));
            mac_file.close();
        }
    }

    if (write) {
        std::ofstream log_file(top_dir/"log.txt");

        log_file << "Crossbar size: " << M << "x" << N << std::endl;
        log_file << "Weight ratio: " << weight_ratio << std::endl;
        log_file << "Input ratio: " << input_ratio << std::endl;
        log_file << "Batch size: " << batch_size << std::endl;
        log_file << "Average execution time: " << total_time/batch_size << " (ms)" << std::endl;
        log_file.close();
    }

    std::cout << "\rSimulations finished with " << total_time/batch_size << " (ms) average execution time" << std::endl;
}

