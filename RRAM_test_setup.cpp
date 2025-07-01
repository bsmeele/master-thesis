#include "nonlinear_crossbar_solver.h"
#include "crossbar_model/linear_crossbar_solver.h"
#include "crossbar_simulator.h"
#include "simulation_settings.h"
#include "memristor.h"
#include "JART_VCM_v1b_var.h"
#include "dummy_memristor.h"
#include "dummy_nonlinear_memristor.h"
#include "VTEAM.h"

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
#include <cmath>

int main(int argc, char* argv[]) {
    int M = 64;
    int N = 64;
    int batch_size = 1;
    bool write = false;
    float weight_ratio = 0.5;
    float input_ratio = 0.5;
    unsigned int seed = time(0);
    bool linear = false;
    bool drift = true;
    int model = 0;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if ((arg == "-m" || arg== "-M" || arg == "-row" || arg == "-rows") && argc > i+1) {
            M = std::stoi(argv[i+1]);
        } else if ((arg == "-n" || arg == "-N" || arg == "-col" || arg == "-cols") && argc > i+1) {
            N = std::stoi(argv[i+1]);
        } else if ((arg == "-b" || arg == "-batch_size") && argc > i+1) {
            batch_size = std::stoi(argv[i+1]);
        } else if ((arg == "-w" || arg == "-write") && argc > i+1) {
            write = (std::atoi(argv[i+1]) == 1) ? true : false;
        } else if ((arg == "-wr" || arg == "-weight_ratio") && argc > i+1) {
            weight_ratio = std::atof(argv[i+1]);
        } else if ((arg == "-ir" || argv[i] == "-input_ratio") && argc > i+1) {
            input_ratio = std::atof(argv[i+1]);
        } else if ((arg == "-s" || arg == "-seed") && argc > i+1) {
            seed = std::stoul(argv[i+1]);
        } else if ((arg == "-l" || arg == "-linear") && argc > i+1) {
            linear = (std::atoi(argv[i+1]) == 1) ? true : false;
        } else if ((arg == "-d" || arg == "-drift") && argc > i+1) {
            drift = (std::atoi(argv[i+1]) == 1) ? true : false;
        } else if ((arg == "-model") && argc > i+1) {
            model = std::stoi(argv[i+1]);
        }
    }

    std::cout << "Random seed: " << seed << std::endl;

    srand(seed);
    // srand((unsigned int) time(0));
    // srand(42);

    // int M = (argc >= 2) ? std::atoi(argv[1]) : 64;
    // int N = (argc >= 3) ? std::atoi(argv[2]) : 64;

    // int batch_size = (argc >= 4) ? std::atoi(argv[3]) : 1;

    // bool write = (argc >= 5 && std::atoi(argv[4]) == 1) ? true : false;

    // float weight_ratio = (argc >= 6) ? std::atof(argv[5]) : 0.5;  // 0 is all HRS, 1 is all LRS
    // float input_ratio = (argc >= 7) ? std::atof(argv[6]) : 0.5;  // 0 is all off, 1 is all on

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
    switch(model) {
        case 0:
            crossbar.Initialize<JART_VCM_v1b_var>();
            break;
        case 1:
            crossbar.Initialize<Dummy>();
            break;
        case 2:
            crossbar.Initialize<DummyNonlinear>();
            break;
        case 3:
            crossbar.Initialize<VTEAM>();
            break;
        default:
            crossbar.Initialize<JART_VCM_v1b_var>();
    }
    
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
        ss << i << " out of " << batch_size << " sims completed (" << std::fixed << std::setprecision(1) << (float) i/batch_size * 100. << "%) ";
        std::cout << "\r" << ss.str() << std::flush;
        
        // Randomly initialize the weigths based on the ratio
        for (int m = 0; m < M; m++) {
            for (int n = 0; n < N; n++) {
                weights[m][n] = ((double) rand())/RAND_MAX < weight_ratio;
            }
        }

        // Randomly initialize the input based on the ratio
        for (int m = 0; m < M; m++) {
            input[m] = ((double) rand())/RAND_MAX < input_ratio;
        }

        std::vector<std::vector<float>> Iout_avg;
        std::vector<float> Iout_MAC;
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        crossbar.Simulate(input, {}, {}, {}, weights, Vwave, dt, Iout_avg, Iout_MAC, linear, drift);
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto execution_time = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
        total_time += execution_time;

        // Write all to files
        if (write) {
            int64_t M64 = static_cast<int64_t>(M);
            int64_t N64 = static_cast<int64_t>(N);

            std::string sim_dir = "sim_" + std::to_string(i);
            std::filesystem::create_directory(top_dir/sim_dir);

            std::ofstream weight_file(top_dir/sim_dir/"weight.bin", std::ios::binary);
            weight_file.write(reinterpret_cast<const char*>(&M64), sizeof(int64_t));
            weight_file.write(reinterpret_cast<const char*>(&N64), sizeof(int64_t));
            for (int m = 0; m < M; m++) {
                for (int n = 0; n < N; n++) {
                    int64_t value = weights[m][n] ? 1 : 0;
                    weight_file.write(reinterpret_cast<const char*>(&value), sizeof(int64_t));
                }
            }
            weight_file.close();

            std::ofstream input_file(top_dir/sim_dir/"input.bin", std::ios::binary);
            input_file.write(reinterpret_cast<const char*>(&M64), sizeof(int64_t));
            for (int m = 0; m < M; m++) {
                int64_t value = input[m] ? 1 : 0;
                input_file.write(reinterpret_cast<const char*>(&value), sizeof(int64_t));
            }
            input_file.close();

            std::ofstream out_file(top_dir/sim_dir/"out.bin", std::ios::binary);
            out_file.write(reinterpret_cast<const char*>(&M64), sizeof(int64_t));
            out_file.write(reinterpret_cast<const char*>(&N64), sizeof(int64_t));
            for (int m = 0; m < M; m++) {
                out_file.write(reinterpret_cast<const char*>(Iout_avg[m].data()), N * sizeof(float));
            }
            out_file.close();
            
            std::ofstream mac_file(top_dir/sim_dir/"out_MAC.bin", std::ios::binary);
            mac_file.write(reinterpret_cast<const char*>(&M64), sizeof(int64_t));
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
        log_file << "Random seed: " << seed << std::endl;
        log_file << "Linear: " << linear << std::endl;
        log_file << "Average execution time: " << (total_time/1000.)/batch_size << " (ms)" << std::endl;

        log_file << std::endl << "Voltage pulse parameters: " << std::endl;
        log_file << "voltage_pulse_width: " << voltage_pulse_width << std::endl;
        log_file << "voltage_pulse_height: " << voltage_pulse_height << std::endl;
        log_file << "voltage_pulse_rise_time: " << voltage_pulse_rise_time << std::endl;
        log_file << "voltage_pulse_fall_time: " << voltage_pulse_fall_time << std::endl;


        log_file << std::endl << "Memristor model paramerters: " << std::endl;
        log_file << crossbar.RRAM[0][0]->GetParams();

        log_file.close();

        std::cout << "\rSaved to: " << top_dir << std::endl;
    }

    std::cout << "\rSimulations finished with " << (total_time/1000.)/batch_size << " (ms) average execution time" << std::endl;
}
