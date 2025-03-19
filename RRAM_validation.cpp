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

std::vector<std::vector<std::vector<float>>> readTensorFromFile(const std::filesystem::path& filePath) {
    std::ifstream file(filePath, std::ios::binary);

    if (!file) {
        throw std::runtime_error("Failed to open file: " + filePath.string());
    }

    // Read dimensions (depth, rows, cols)
    int64_t depth, rows, cols;
    file.read(reinterpret_cast<char*>(&depth), sizeof(int64_t));
    file.read(reinterpret_cast<char*>(&rows), sizeof(int64_t));
    file.read(reinterpret_cast<char*>(&cols), sizeof(int64_t));

    // Create a 3D tensor
    std::vector<std::vector<std::vector<float>>> tensor(depth, std::vector<std::vector<float>>(rows, std::vector<float>(cols)));

    // Read the data into the tensor
    for (int64_t i = 0; i < depth; ++i) {
        for (int64_t j = 0; j < rows; ++j) {
            file.read(reinterpret_cast<char*>(tensor[i][j].data()), cols * sizeof(float));
        }
    }

    return tensor;
}

std::vector<std::vector<float>> readFloatMatrixFromFile(const std::filesystem::path& filePath) {
    std::ifstream file(filePath, std::ios::binary);

    if (!file) {
        throw std::runtime_error("Failed to open file: " + filePath.string());
    }

    // Read dimensions (rows, cols)
    int64_t rows, cols;
    file.read(reinterpret_cast<char*>(&rows), sizeof(int64_t));
    file.read(reinterpret_cast<char*>(&cols), sizeof(int64_t));

    // Create a 2D matrix
    std::vector<std::vector<float>> matrix(rows, std::vector<float>(cols));

    // Read the data into the matrix
    for (int64_t i = 0; i < rows; ++i) {
        file.read(reinterpret_cast<char*>(matrix[i].data()), cols * sizeof(float));
    }

    return matrix;
}

std::vector<std::vector<double>> readMac(const std::filesystem::path& filePath) {
    std::ifstream file(filePath, std::ios::binary);

    if (!file) {
        throw std::runtime_error("Failed to open file: " + filePath.string());
    }

    int64_t shape[2];
    file.read(reinterpret_cast<char*>(shape), 2 * sizeof(int64_t)); 

    int64_t rows = shape[0];
    int64_t cols = shape[1];

    // Read the matrix data
    // std::vector<double> data(rows * cols);
    // file.read(reinterpret_cast<char*>(data.data()), rows * cols * sizeof(double));

    std::vector<std::vector<double>> data(rows, std::vector<double>(cols));
    // Read the data into the matrix
    for (int64_t i = 0; i < rows; ++i) {
        file.read(reinterpret_cast<char*>(data[i].data()), cols * sizeof(double));
    }

    if (!file) {
        throw std::runtime_error("Error reading matrix data.");
    }

    return data;
}

std::vector<std::vector<std::vector<double>>> readMem(const std::filesystem::path& filePath) {
    std::ifstream file(filePath, std::ios::binary);

    if (!file) {
        throw std::runtime_error("Failed to open file: " + filePath.string());
    }

    int64_t shape[3];
    file.read(reinterpret_cast<char*>(shape), 3 * sizeof(int64_t));

    int64_t depth = static_cast<int>(shape[0]);
    int64_t rows = static_cast<int>(shape[1]);
    int64_t cols = static_cast<int>(shape[2]);

    std::vector<std::vector<std::vector<double>>> data(depth, std::vector<std::vector<double>>(rows, std::vector<double>(cols)));
    // Read the data into the matrix
    for (int64_t i = 0; i < depth; ++i) {
        for (int64_t j = 0; j < rows; j++) {
                file.read(reinterpret_cast<char*>(data[i][j].data()), cols * sizeof(double));
        }
    }

    if (!file) {
        throw std::runtime_error("Error reading matrix data.");
    }

    return data;
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

    // 1 ohm resistors
    // Ndiscmin of 0.0001, default parameters otherwise
    // Transistor gate is connected to source line, thus memristor will be connected if source line is 1
    CrossbarSimulator crossbar = CrossbarSimulator(M, N);
    // crossbar.Rswl1 = 0.;
    // crossbar.Rswl2 = 0.;
    // crossbar.Rsbl1 = 0.;
    // crossbar.Rsbl2 = 0.;
    // crossbar.Rwl = 0.;
    // // crossbar.Rbl = 0.;
    // for (int i = 0; i < crossbar.RRAM.size(); i++) {
    //     for (int j = 0; j < crossbar.RRAM[i].size(); j++) {
    //         crossbar.RRAM[i][j].Ndiscmin = 0.0001;
    //         crossbar.RRAM[i][j].Ninit = 0.0001;
    //     }
    // }

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
    // for (int i = 0; i < M; i++) {
    //     std::vector<bool> row;
    //     for (int j = 0; j < N; j++) {
    //         if (std::rand()%2 == 0) {
    //             row.push_back(true);
    //         } else {
    //             row.push_back(false);
    //         }
    //     }
    //     weights.push_back(row);
    // }

    // std::vector<bool> Vapp;
    // for (int i = 0; i < M; i++) {
    //     if (std::rand()%2 == 0) {
    //         Vapp.push_back(true);
    //     } else {
    //         Vapp.push_back(false);
    //         crossbar.access_transistors[i] = std::vector<bool>(N, false);
    //     }
    // }

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

        // 0.1 V read voltage
        // width of 50 us
        // 5 us rise and fall time
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
        // std::vector<std::vector<float>> Iout_avg = Iwave[(voltage_pulse_width/simulation_time_step)/2.];

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

    // for (int i = 0; i < 32; i++) {
    //     std::cout << output_data_MAC[0][i] << " ";
    // }
    // std::cout << std::endl;

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

    // for (int i = 0; i < output_data.size(); i++) {
    //     for (int m = 0; m < output_data[i].size(); m++) {
    //         for (int n = 0; n < output_data[i][m].size(); n++) {
    //             outfile << output_data[i][m][n] << " ";
    //         }
    //         outfile << std::endl;
    //     }
    //     outfile << std::endl;
    // }

    outfile.close();

    outfile = std::ofstream(top_dir/"out_mac.bin", std::ios::binary);

    if (!outfile) {
        std::cout << "Failed to open file: " << top_dir/"out_mac.bin" << std::endl;
        return 1;
    }

    rows = output_data_MAC.size();
    cols = output_data_MAC.empty() ? 0 : output_data_MAC[0].size();

    outfile.write(reinterpret_cast<const char*>(&rows), sizeof(int64_t));
    outfile.write(reinterpret_cast<const char*>(&cols), sizeof(int64_t));

    for (int64_t i = 0; i < rows; ++i) {
        outfile.write(reinterpret_cast<const char*>(output_data_MAC[i].data()), cols * sizeof(float));
    }

    // for (int i = 0; i < output_data_MAC.size(); i++) {
    //     for (int j = 0; j < output_data_MAC[i].size(); j++) {
    //         outfile << output_data_MAC[i][j] << " ";
    //     }
    //     outfile << std::endl;
    // }

    outfile.close();
    
    // std::vector<std::vector<double>> mac_data;
    // try { mac_data = readMac(mac_path); }
    // catch (const std::exception& e) {
    //     // std::cout << e.what() << std::endl;
    //     std::cout << "No validation data" << std::endl;
    //     return 0;
    // }

    // std::vector<std::vector<std::vector<double>>> mem_data;
    // try { mem_data = readMem(mem_path); }
    // catch (const std::exception& e) {
    //     // std::cout << e.what() << std::endl;
    //     std::cout << "No validation data" << std::endl;
    //     return 0;
    // }

    // assert(mac_data.size() == output_data_MAC.size());
    // assert(mac_data[0].size() == output_data_MAC[0].size());
    // assert(mem_data.size() == output_data.size());
    // assert(mem_data[0].size() == output_data[0].size());
    // assert(mem_data[0][0].size() == output_data[0][0].size());

    // for (int i = 0; i < input_data[0].size(); i++) {
    //     std::cout << input_data[0][i] << " ";
    // }
    // std::cout << std::endl << std::endl;

    // for (int i = 0; i < weight_data.size(); i++) {
    //     for (int j = 0; j < weight_data[0].size(); j++) {
    //         std::cout << weight_data[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    // float err = 0.;
    // for (int i = 0; i < mac_data[0].size(); i++) {
    //     // std::cout << mac_data[0][i] * 1e-6 - output_data_MAC[0][i] << " ";
    //     err += fabs(mac_data[0][i] * 1e-6 - output_data_MAC[0][i]);
    // }
    // // std::cout << std::endl;
    // err = err / mac_data[0].size();
    // std::cout << "Average mac error: " << err << std::endl;

    // err = 0.;
    // for (int i = 0; i < mem_data[0].size(); i++) {
    //     for (int j = 0; j < mem_data[0][i].size(); j++) {
    //         // std::cout << mem_data[0][i][j] * 1e-6 - output_data[0][i][j] << " ";
    //         err += fabs(mem_data[0][i][j] * 1e-6 - output_data[0][i][j]);
    //     }
    //     // std::cout << std::endl;
    // }
    // err = err / (mem_data[0].size() * mem_data[0][0].size());
    // std::cout << "Average individual error: " << err << std::endl;

    // for (int i = 0; i < mem_data[0].size(); i++) {
    //     for (int j = 0; j < mem_data[0][i].size(); j++) {
    //         std::cout << mem_data[0][i][j] * 1e-6 << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    // for (int i = 0; i < mem_data[0].size(); i++) {
    //     for (int j = 0; j < mem_data[0][i].size(); j++) {
    //         std::cout << output_data[0][i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // std::vector<float> Iout;
    // for (int n = 0; n < N; n++) {
    //     float Iavg = 0;
    //     for (int i = 0; i < 40; i++) {
    //         for (int m = 0; m < M; m++) {
    //             Iavg += Iwave[i][m][n];
    //         }
    //     }
    //     Iavg /= 40.;
    //     Iout.push_back(Iavg);
    // }

    // if (print) {
    //     std::cout << "Vapplied (boolean):" << std::endl;
    //     for (int i = 0; i < Vapp.size(); i++) {
    //         std::cout << Vapp[i] << std::endl;
    //     }
    //     std::cout << std::endl;

    //     std::cout << "Weights (boolean):" << std::endl;
    //     for (int i = 0; i < weights.size(); i++) {
    //         for (int j = 0; j < weights[0].size(); j++) {
    //             std::cout << weights[i][j] << " ";
    //         }
    //         std::cout << std::endl;
    //     }
    //     std::cout << std::endl;

    //     std::cout << "Iout:" << std::endl;
    //     for (int i = 0; i < Iout.size(); i++) {
    //         std::cout << Iout[i] << std::endl;
    //     }
    // }

    // std::ofstream outfile("out.txt");

    // if (!outfile) {
    //     std::cout << "No out file" << std::endl;
    //     return 1;
    // }

    // for (int i = 0; i < Vapp.size(); i++) {
    //     outfile << Vapp[i] << std::endl;
    // }
    // outfile << std::endl;

    // for (int i = 0; i < weights.size(); i++) {
    //     for (int j = 0; j < weights[0].size(); j++) {
    //         outfile << weights[i][j] << " ";
    //     }
    //     outfile << std::endl;
    // }
    // outfile << std::endl;

    // for (int i = 0; i < M; i++) {
    //     for (int j = 0; j < N; j++) {
    //         outfile << Iwave[25][i][j] << " ";
    //     }
    //     outfile << std::endl;
    // }
    // outfile << std::endl;

    // outfile.close();
}
