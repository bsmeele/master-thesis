#include "memristor_model/JART_VCM_v1b_var.h"
#include "crossbar_model/cam.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <array>

// The crossbar solver assumes all devices to be linear. However, the memristor model is nonlinear
// This means some adaptations are required to combine the two
// After promting chatGPT, I've come up with some options:
//   Use an outer iterative solver to make a guess at the nodal voltages, and thus the associated memristor conductances
//   Assume linearity, perhaps by making a guess of the nodal votlages and assuming linearity around that point (doesn't seem very promising)
//   Use look up tables in some way by making a matrix with some select votlage-state pairs and using interpolation for in between values
// For now, an outer iterative solver seems like the best solution
// Several papers have mentioned the Newton-Raphson method for this. However, this requires a Jacobian matrix (a matrix containing all partial derivatives)
//   and I have no idea how to go about that
// For the sake of getting results in a reasonable amount of time, I think a fixed point solver will have to be sufficient
// The fixed point solver will work as follows:
//   Make some initial guess of the nodal votlages (using the results of a previous state might speed things up)
//   Determine the conductance matrix based on the guess
//   Use the conductance matrix to calculated the associated nodal voltages
//   Check if the calculated nodal voltages are equal (or close) to the guess
//   If not, the calculated voltages becomes the new guess.
//   Repeat untill convergence

// Initial fixed point implementation:
//   Very slow (5 seconds for 3x3), it reaches 1e-2 pretty fast, but takes exponentially longer to reach 1e-3 and lower
//   Adding a relaxation factor (V = w * Vnew + (1-w) * Vold) of 1.3 seemed to speed things up somewhat (3.8s for 3x3) but the main issue persists
//   Using a relaxation factor of 1.4 makes it so a NAN is encountered in a memristor, suggesting using relaxation like this might be unstable (or I need to look into why it's producing NANs)
//     Making the relaxation factor too high seems to make the norm worse
//   Similarly, using a dynamic relaxation factor (starting at 1.1, adding 0.0001 every iteration) slightly increases the speed (3.6s for 3x3), but making the dynamic increase too high makes it unstable

// I no longer think fixed point is sufficient. So here is a primer on the Newton-Raphson method:
//   This method is an iterative method to solve F(V) = 0
//   It works similarly to fixed point in that it starts with a guess and updates that guess, but it updates the guess as such: Vnew = Vold - J^(-1)(Vold) * F(Vold)
//   F(V) in this case is the residual function: Vnew - Vold
//   J(V) is the Jacobian matrix of F(V), the matrix containing every partial derivative of F(V)
//   J(V) is subsequently also the main difficulty of this method. At least for this problem, the Jacobain is hard to define
//   One option is approximating it with finite differences: J(V) = (F(V + dV) - F(V)) / dV
//   Another option is iteratively update it using Broyden's method: Jnew = Jold + ((Fnew - Fold - Jold * dV) * dv^T) / ||dV^2||
//   For now I think I'll be using finite differences
//   Some alterations:
//     Solve J(V) * dV = -F(V) with an iterative solver, and update V with Vnew = Vold + dV
//     Use a dampening factor to avoid overshooting: Vnew = Vold + a * dV

Eigen::VectorXf newton_raphson_solve(
    std::vector<std::vector<JART_VCM_v1b_var>> RRAM,
    Eigen::VectorXf& Vguess, Eigen::SparseMatrix<float>& G_ABCD,
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
    const float Rswl1, const float Rswl2, const float Rsbl1, const float Rsbl2,
    const float Rwl, const float Rbl,
    const bool print = false
) {
    if (RRAM.size() == 0) { return Eigen::VectorXf(0); }
    int M = RRAM.size();
    int N = RRAM[0].size();

    Eigen::MatrixXf G(M, N);

    int it = 0;
    while(true) {// Determine G
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                float v = Vguess(i*N + j) - Vguess(i*N + j + M*N);
                if (v == 0) { v += 1e-6; }
                G(i, j) = RRAM[i][j].apply_voltage(v, 0) / v;
            }
        }

        // Calculate V
        Eigen::VectorXf Vout = solve_cam(G, Vguess, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);

        Eigen::VectorXf Fv = Vout - Vguess;

        std::cout << "Norm: " << Fv.norm() << std::endl;
        if (Fv.norm() < 1e-6) {
            if (print) {
                std::cout << "V:\n" << Vout << std::endl << std::endl;
                std::cout << "G:\n" << G << std::endl << std::endl;
                std::cout << "R:\n" << G.array().inverse() << std::endl << std::endl;
                std::cout << "Solved in " << it << " iterations" << std::endl;
            }
            return Vout;
        }

        // Calculate Jacobian
        Eigen::MatrixXf J = Eigen::MatrixXf::Zero(2*M*N, 2*M*N);
        float delta = 1e-6;
        for (int col = 0; col < 2*M*N; col++) {
            std::cout << "Col: " << col << std::endl;

            Vguess(col) += delta;
            // std::cout << Vguess << std::endl << std::endl;

            for (int i = 0; i < M; i++) {
                for (int j = 0; j < N; j++) {
                    float v = Vguess(i*N + j) - Vguess(i*N + j + M*N);
                    std::cout << "(i, j): " << i << " " << j << std::endl;
                    if (v == 0) { v += 1e-6; }
                    G(i, j) = RRAM[i][j].apply_voltage(v, 0) / v;
                }
            }

            Eigen::VectorXf Fj = solve_cam(G, Vguess, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl) - Vguess;
            J.col(col) = (Fj - Fv) / delta;

            Vguess(col) -= delta;
        }

        // Calculate dV
        Eigen::VectorXf dV = J.partialPivLu().solve(-Fv);

        // Update guess
        Vguess = Vguess + dV;

        it++;
    }
}

Eigen::VectorXf fixedpoint_solve(
    std::vector<std::vector<JART_VCM_v1b_var>> RRAM,
    Eigen::VectorXf& Vguess, Eigen::SparseMatrix<float>& G_ABCD,
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
    const float Rswl1, const float Rswl2, const float Rsbl1, const float Rsbl2,
    const float Rwl, const float Rbl,
    const bool print = false
) {
    if (RRAM.size() == 0) { return Eigen::VectorXf(0); }
    int M = RRAM.size();
    int N = RRAM[0].size();

    Eigen::MatrixXf G(M, N);

    // float w = 1.1;

    int it = 0;
    while (true) {
        // std::cout << V << std::endl << std::endl;

        // Determine G
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                float v = Vguess(i*N + j) - Vguess(i*N + j + M*N);
                if (v == 0) { v += 1e-6; }
                G(i, j) = RRAM[i][j].apply_voltage(v, 0) / v;
            }
        }
        
        // Calculate V
        Eigen::VectorXf Vout = solve_cam(G, Vguess, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);

        if (print) { std::cout << "Norm: " << (Vout - Vguess).norm() << std::endl; }
        // Exit condition
        if ((Vout - Vguess).norm() < 1e-3) {
            if (print) {
                std::cout << "V:\n" << Vout << std::endl << std::endl;
                std::cout << "G:\n" << G << std::endl << std::endl;
                std::cout << "R:\n" << G.array().inverse() << std::endl << std::endl;
                std::cout << "Solved in " << it << " iterations" << std::endl;
            }
            return Vout;
        }

        Vguess = Vout;

        // V = w * V + (1 - w) * V_guess;

        // w += 0.0001;

        it++;
    }
}

int main(int argc, char* argv[]) {
    srand((unsigned int) time(0));

    int M = (argc >= 2) ? std::atoi(argv[1]) : 16;
    int N = (argc >= 3) ? std::atoi(argv[2]) : 16;

    int runs = (argc >= 4) ? std::atoi(argv[3]) : 1;

    bool print = (argc >= 5 && std::atoi(argv[4]) == 1) ? true : false;
    
    Eigen::VectorXf V = Eigen::VectorXf::Zero(2*M*N);

    float Rswl1 = 3.;
    float Rswl2 = INFINITY;
    float Rsbl1 = INFINITY;
    float Rsbl2 = 5.;

    float Rwl = 3.;
    float Rbl = 2.;

    Eigen::SparseMatrix<float> G_ABCD = partially_precompute_G_ABCD(M, N, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);

    std::vector<std::vector<JART_VCM_v1b_var>> RRAM;
    for (int i = 0; i < M; i++) {
        std::vector<JART_VCM_v1b_var> row;
        for (int j = 0; j < N; j++) {
            row.push_back(JART_VCM_v1b_var());
        }
        RRAM.push_back(row);
    }

    float Vdd = 5.;
    Eigen::VectorXf Vappwl1 = Eigen::VectorXf::Random(M);
    Vappwl1 = (Vappwl1.array() > 0.5).select(Eigen::VectorXf::Constant(M, Vdd), Eigen::VectorXf::Zero(M));
    Vappwl1(0) = 5;
    Vappwl1(1) = 7;
    Vappwl1(2) = 9;

    Eigen::VectorXf Vappwl2 = Eigen::VectorXf::Zero(M);
    Eigen::VectorXf Vappbl1 = Eigen::VectorXf::Zero(M);
    Eigen::VectorXf Vappbl2 = Eigen::VectorXf::Zero(M);

    if (print) {
        std::cout << "Vappwl1:\n" << Vappwl1 << std::endl << std::endl;
    }

    // {
    // auto start_time = std::chrono::high_resolution_clock::now();

    // Eigen::VectorXf Vout = fixedpoint_solve(RRAM, V, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, print);

    // auto end_time = std::chrono::high_resolution_clock::now();

    // auto execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    // std::cout << "Execution time: " << execution_time << " (ms)" << std::endl;

    // if (print) {
    //     Eigen::MatrixXf G(M, N);
    //     for (int i = 0; i < M; i++) {
    //         for (int j = 0; j < N; j++) {
    //             float v = V(i*N + j) - V(i*N + j + M*N);
    //             if (v == 0) { v += 1e-6; }
    //             G(i, j) = RRAM[i][j].apply_voltage(v, 0) / v;
    //         }
    //     }

    //     std::vector<float> Iout;
    //     for (int j = 0; j < N; j++) {
    //         float Ioutj = 0;
    //         for (int i = 0; i < M; i++) {
    //             Ioutj += (Vout(i*N + j) - Vout(i*N + j + M*N)) * G(i,j);
    //         }
    //         Iout.push_back(Ioutj);
    //     }

    //     std::cout << "Iout:" << std::endl;
    //     for (int j = 0; j < N; j++) {
    //         std::cout << Iout[j] << std::endl;
    //     }
    //     std::cout << std::endl;
    // }
    // }

    {
    auto start_time = std::chrono::high_resolution_clock::now();

    Eigen::VectorXf Vout = newton_raphson_solve(RRAM, V, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, print);

    auto end_time = std::chrono::high_resolution_clock::now();

    auto execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "Execution time: " << execution_time << " (ms)" << std::endl;

    if (print) {
        Eigen::MatrixXf G(M, N);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                float v = V(i*N + j) - V(i*N + j + M*N);
                if (v == 0) { v += 1e-6; }
                G(i, j) = RRAM[i][j].apply_voltage(v, 0) / v;
            }
        }

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
    }
}
