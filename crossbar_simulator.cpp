#include "crossbar_simulator.h"

#include "nonlinear_crossbar_solver.h"
#include "crossbar_model/linear_crossbar_solver.h"
#include "simulation_settings.h"

#include <iostream>

void CrossbarSimulator::SetRRAM(std::vector<std::vector<bool>> weights) {
    assert(weights.size() == RRAM.size());
    assert(weights[0].size() == RRAM[0].size());

    for (int i = 0; i < weights.size(); i++) {
        for (int j = 0; j < weights[0].size(); j++) {
            RRAM[i][j]->SetWeight(weights[i][j]);
        }
    }
}

void CrossbarSimulator::SetAccessTransistors(std::vector<bool> gate_lines) {
    assert(gate_lines.size() == RRAM.size());
    for (int m = 0; m < M; m++) {
        if (gate_lines[m]) {
            access_transistors[m] = std::vector<bool>(N, true);
        } else {
            access_transistors[m] = std::vector<bool>(N, false);
        }
    }
}

Eigen::VectorXf CrossbarSimulator::LinearSolve(
        Eigen::VectorXf Vguess,
        const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
        const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2
) {
    Eigen::VectorXf E = ComputeE(M, N, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);

    Eigen::MatrixXf G(M, N);
    if (simulation_num_threads > 0) {
        std::vector<std::future<void>> futures;
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                if (access_transistors[i][j]) {
                    futures.emplace_back(pool.enqueue([&, i, j]() {
                        float v = Vguess(i*N + j) - Vguess(i*N + j + M*N);
                        G(i, j) = (float) 1./RRAM[i][j]->GetResistance(v);
                    }));
                } else {
                    G(i, j) = 0;
                }
            }
        }

        for (auto& fut : futures) {
            fut.get();
        } 
    } else {
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                if (access_transistors[i][j]) {
                    float v = Vguess(i*N + j) - Vguess(i*N + j + M*N);
                    G(i, j) = (float) 1./RRAM[i][j]->GetResistance(v);
                } else {
                    G(i, j) = 0;
                }
            }
        }
    }

    return SolveCam(G, Vguess, partial_G_ABCD, E, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, linear_solver);
}

Eigen::VectorXf CrossbarSimulator::NonlinearSolve(
    Eigen::VectorXf Vguess,
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
    std::string method
) {
    Eigen::VectorXf E = ComputeE(M, N, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);

    if (method == "fixed-point") {
        return FixedpointSolve(RRAM, access_transistors, Vguess, partial_G_ABCD, E, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, linear_solver, pool);
    } else if (method == "NewtonRaphson") {
        return NewtonRaphsonSolve(RRAM, access_transistors, Vguess, partial_G_ABCD, E, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, linear_solver);
    } else if (method == "Broyden") {
        return BroydenSolve(RRAM, access_transistors, Vguess, partial_G_ABCD, E, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, linear_solver);
    } else if (method == "BroydenInv") {
        return BroydenInvSolve(RRAM, access_transistors, Vguess, partial_G_ABCD, E, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, linear_solver);
    } else {
        return FixedpointSolve(RRAM, access_transistors, Vguess, partial_G_ABCD, E, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, linear_solver, pool);
    }
}

std::vector<std::vector<float>> CrossbarSimulator::ApplyVoltage(
    Eigen::VectorXf Vguess,
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
    float dt, std::string method
) {
    Eigen::VectorXf Vout = NonlinearSolve(Vguess, Vappwl1, Vappwl2, Vappbl1, Vappbl2, method);

    std::vector<std::vector<float>> Iout(M, std::vector<float>(N));

    if (simulation_num_threads > 0) {
        std::vector<std::future<void>> futures;
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                if (access_transistors[i][j]) {
                    futures.emplace_back(pool.enqueue([&, i, j]() {
                        float v = Vout(i*N + j) - Vout(i*N + j + M*N);
                        Iout[i][j] = RRAM[i][j]->ApplyVoltage(v, dt);
                    }));
                } else {
                    Iout[i][j] = 0.;
                }
            }
        }

        for (auto& fut : futures) {
            fut.get();
        } 
    } else {
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                if (access_transistors[i][j]) {
                    float v = Vout(i*N + j) - Vout(i*N + j + M*N);
                    Iout[i][j] = RRAM[i][j]->ApplyVoltage(v, dt);
                } else {
                    Iout[i][j] = 0.;
                }
            }
        }
    }

    return Iout;
}

std::vector<float> CrossbarSimulator::CalculateIout(Eigen::VectorXf Vout) {
    std::vector<float> Iout;
    for (int j = 0; j < N; j++) {
        float Ioutj = 0.;
        for (int i = 0; i < M; i++) {
            float v = Vout(i*N + j) - Vout(i*N + j + M*N);
            Ioutj += RRAM[i][j]->ApplyVoltage(v, 0);
        }
        Iout.push_back(Ioutj);
    }
    return Iout;
}

void CrossbarSimulator::SetCrossbarParameters(float Rswl1, float Rswl2, float Rsbl1, float Rsbl2, float Rwl, float Rbl) {
    this->Rswl1 = Rswl1;
    this->Rswl2 = Rswl2;
    this->Rsbl1 = Rsbl1;
    this->Rsbl2 = Rsbl2;
    
    this->Rwl = Rwl;
    this->Rbl = Rbl;
    
    UpdatePrecomputeG_ABCD();
}

void CrossbarSimulator::UpdatePrecomputeG_ABCD() {
    partial_G_ABCD = PartiallyPrecomputeG_ABCD(M, N, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);
}
