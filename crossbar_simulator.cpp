#include "crossbar_simulator.h"

#include "nonlinear_crossbar_solver.h"
#include "crossbar_model/linear_crossbar_solver.h"

#include <iostream>

void CrossbarSimulator::SetRRAM(std::vector<std::vector<bool>> weights) {
    assert(weights.size() == RRAM.size());
    assert(weights[0].size() == RRAM[0].size());

    for (int i = 0; i < weights.size(); i++) {
        for (int j = 0; j < weights[0].size(); j++) {
            if (weights[i][j]) { RRAM[i][j].Nreal = RRAM[i][j].Ndiscmax; }
            else { RRAM[i][j].Nreal = RRAM[i][j].Ndiscmin; }
            RRAM[i][j].Treal = RRAM[i][j].T0;
        }
    }
}

Eigen::VectorXf CrossbarSimulator::NonlinearSolve(
    Eigen::VectorXf Vguess,
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
    std::string method
) {
    Eigen::VectorXf E = ComputeE(M, N, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);

    if (method == "fixed-point") {
        return FixedpointSolve(RRAM, access_transistors, Vguess, partial_G_ABCD, E, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, linear_solver);
    } else if (method == "NewtonRaphson") {
        return NewtonRaphsonSolve(RRAM, access_transistors, Vguess, partial_G_ABCD, E, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, linear_solver);
    } else if (method == "Broyden") {
        return BroydenSolve(RRAM, access_transistors, Vguess, partial_G_ABCD, E, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, linear_solver);
    } else if (method == "BroydenInv") {
        return BroydenInvSolve(RRAM, access_transistors, Vguess, partial_G_ABCD, E, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, linear_solver);
    } else {
        return FixedpointSolve(RRAM, access_transistors, Vguess, partial_G_ABCD, E, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, linear_solver);
    }
}

std::vector<std::vector<float>> CrossbarSimulator::ApplyVoltage(
    Eigen::VectorXf Vguess,
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
    float dt, std::string method
) {
    Eigen::VectorXf Vout = NonlinearSolve(Vguess, Vappwl1, Vappwl2, Vappbl1, Vappbl2, method);

    std::vector<std::vector<float>> Iout;
    for (int i = 0; i < M; i++) {
        std::vector<float> row;
        for (int j = 0; j < N; j++) {
            if (!access_transistors[i][j]) {
                row.push_back(0.);
                continue;
            }
            float v = Vout(i*N + j) - Vout(i*N + j + M*N);
            float I = RRAM[i][j].ApplyVoltage(v, dt);
            row.push_back(I);
        }
        Iout.push_back(row);
    }

    return Iout;
}

std::vector<float> CrossbarSimulator::CalculateIout(Eigen::VectorXf Vout) {
    std::vector<float> Iout;
    for (int i = 0; i < M; i++) {
        float Ioutj;
        for (int j = 0; j < N; j++) {
            float v = Vout(i*N + j) - Vout(i*N + j + M*N);
            Ioutj += v / RRAM[i][j].GetResistance(v);
        }
        Iout.push_back(Ioutj);
    }
    return Iout;
}
