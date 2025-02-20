#include "crossbar_simulator.h"

#include "nonlinear_crossbar_solver.h"
#include "crossbar_model/linear_crossbar_solver.h"

void CrossbarSimulator::SetRRAM(std::vector<std::vector<bool>> weights) {
    assert(weights.size() == RRAM.size());
    assert(weights[0].size() == RRAM[0].size());

    for (int i = 0; i < weights.size(); i++) {
        for (int j = 0; j < weights[0].size(); j++) {
            if (weights[i][j]) { RRAM[i][j].Nreal = RRAM[i][j].Ndiscmax; }
            else { RRAM[i][j].Nreal = RRAM[i][j].Ndiscmin; }
        }
    }
}

std::vector<float> CrossbarSimulator::NonlinearSolve(
    Eigen::VectorXf Vguess,
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
    std::string method
) {
    Eigen::VectorXf Vout;
    if (method == "fixed-point") {
        Vout = FixedpointSolve(RRAM, access_transistors, Vguess, partial_G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);
    } else if (method == "NewtonRaphson") {
        Vout = NewtonRaphsonSolve(RRAM, access_transistors, Vguess, partial_G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);
    } else if (method == "Broyden") {
        Vout = BroydenSolve(RRAM, access_transistors, Vguess, partial_G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);
    } else if (method == "BroydenInv") {
        Vout = BroydenInvSolve(RRAM, access_transistors, Vguess, partial_G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);
    } else {
        Vout = FixedpointSolve(RRAM, access_transistors, Vguess, partial_G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);
    }

    std::vector<float> Iout;
    for (int j = 0; j < N; j++) {
        float Ioutj = 0;
        for (int i = 0; i < M; i++) {
            if (access_transistors[i][j]) {
                float v = Vout(i*N + j) - Vout(i*N + j + M*N);
                Ioutj += v / RRAM[i][j].GetResistance(v);
            }
        }
        Iout.push_back(Ioutj);
    }

    return Iout;
}

std::vector<float> CrossbarSimulator::ApplyVoltage(
    Eigen::VectorXf Vguess,
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
    float dt, std::string method
) {
    Eigen::VectorXf Vout;
    if (method == "fixed-point") {
        Vout = FixedpointSolve(RRAM, access_transistors, Vguess, partial_G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);
    } else if (method == "NewtonRaphson") {
        Vout = NewtonRaphsonSolve(RRAM, access_transistors, Vguess, partial_G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);
    } else if (method == "Broyden") {
        Vout = BroydenSolve(RRAM, access_transistors, Vguess, partial_G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);
    } else if (method == "BroydenInv") {
        Vout = BroydenInvSolve(RRAM, access_transistors, Vguess, partial_G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);
    } else {
        Vout = FixedpointSolve(RRAM, access_transistors, Vguess, partial_G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);
    }

    std::vector<float> Iout;
    for (int i = 0; i < M; i++) {
        float Ioutj = 0;
        for (int j = 0; j < N; j++) {
            float v = Vout(i*N + j) - Vout(i*N + j + M*N);
            Ioutj += RRAM[i][j].ApplyVoltage(v, dt);
        }
        Iout.push_back(Ioutj);
    }

    return Iout;
}
