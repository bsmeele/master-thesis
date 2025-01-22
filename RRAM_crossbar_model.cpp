#include "memristor_model/JART_VCM_v1b_var.h"
#include "crossbar_model/cam.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <array>
#include <limits>
// #include <cmath>

float norm = 0;

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

// Current issues (resolved):
//   The Jacobian seems way too high, resulting in updates that are way too high, resulting in NANs from the memristors (which might be its own problem, but unrelated to what this solver should do)
//   The values of the Jacobian are currently around 1e3/1e4, while Fj - Fv is around 1e-3/1e-2 for a delta of 1e-6
//   The calculation method seems to be correct, and with the delta and the Fj-Fv values, the high Jacobian values do somewhat make sense
//   One option is that the linear solver is sensitive/inaccurate to small changes
//   Or the memristor conductivity is sensitive/inaccurate to small changes
//   My current best guess is with the way I'm calculating the conductivity of the memristors. It is currently calculated as I/V which means it's very unstable around V=0
//   Even with a guard for v=0, it is still very unstable in that region
//   After further studying the code, this might not be the only issue.
//   For now I have set the initial guess to the crossbar voltages, meaning that the solver won't start around the 0 instability (this seems to be a better guess anyway)
//   This makes it so the solver is able to go through a few iterations before it hits a NAN. However, the norm if F(V) is only increasing suggesting further issues
//   So apperantly adjusting the delta to 1e-3 for the Jacobian calculation fixes everything and it works now

// The Newton-Raphson method is implemented but has some issues:
//   It is a lot slower than fixed point
//   It is a little less accurate than fixed point
// These could be all due to the method of Jacobian calculation
// It uses the finite difference method with a delta of 1e-3
// Finite distances is obviously very slow since the Jacobian has many elements
// The delta could also be too large leading to the loss of accuracy for the algorithm, but smaller deltas seem to be unstable
// One possible option to reduce the execution time is to calculate the Jacobian iteratively with something like the Broyden's Method
// Since the Jacobian is calculated with finite differences and a somewhat large delta, it might not even loose that much accuracy
// Turns out Broyden is faster but also somewhat slow, it is also slightly more accurate
// Using an initial Jacobian equal to the identity matrix seems to slightly increase accuracy and increases performance over partial differences
// My best guess as to why it's still somewhat slow is the dV solve method.
// However, Broyden's method should also work with the inverse of the Jacobian, which would resolve this
// I've tested this inverse method and it seems to be faster but lose significant accuracy. Additionally, it does not work (yet) for 16x16 and up

// Currently I have implemented the fixed point, Newton-Raphson, Broyden, and inverse Broyden methods
// Broyden's method seems to be an all round improvement over Newton-Raphson
// It is faster and more accurate but is still not able to simulate past 32x32 due to execution time
// Inverse Broyden's method could solve this, but as of now does not work for 16x16 and up
// Even still, inverse Broyden loses significant accuracy over Broyden's method
// However, all this (in my opinion) does not really matter as fixed point is close to Broyden's method in terms of accuracy and is way faster
// For now I see no reason not to use the fixed point method
// For all methods, generally the norm approaches some minimum and keeps oscillating around that or might degrade some
// The best improvement for now is probably to detect these oscilation so the solver can exit early

// Using any optimiation flags seems to break the code with vdd of 5 (works without optimization). Works for vdd of 1.5
// The NANs seem to come from the crossbar solver
// Compiler optimization also seems to reduce the accuracy

// Actually, NANs also occur without optimization
// I've tested the intput for which it returns NAN (G matrix and Vappwl1) on the crossbar code and it does not return NAN
// Additionally, I've compared it with the associated Vguess and the norm is 0.0012673, which would be correspond to a slight norm reduction
// Ussing the Vguess as the initial guess in the crossbar code also doesn't change anything.
// Thus something about calling the function from the RRAM code breaks it
// Besides fixing this bug. It might be usefull to think about what should happen when the crossbar solver is given an impossible configuration (altough I can't think of any impossible configurations)

// Apperantly using the conjugate gradient solver without guess in the crossbar solver fixes this bug. No idea why. Additionally, the loss of the guess option increases executiont time
// Actually, it seems to fix it for 16x66. For 32x32 it still exists
// Additionally, it now very rarely reduces NAN for 16x16 because sometimes Vapplies is over thousands of volt (positive or negative). Meaning the solver (or model) goes wrong somewhere
// After 1000 runs of 16x16, the inverse broyden method (haven't tested others) does not seem to exhibit this solver bug (the one in the line above), suggesting it is a fault in the fixed point method

// There are two variations: the 'good' and 'bad' Broyden's methods
// The 'bad' method seems way faster, thus has been used for results
// Measurements are from 1 run without optimization
// Often not able to reach criterion
// For 3x3: 2 iterations, 7 ms
// For 8x8: 1.49084e-05 norm, 306 ms
// For 16x16: 1.72063e-05 norm, 1402 ms
// For 32x32: 0.00140215 norm, 13145 ms
// For 64x64: 0.00852754 norm, 158731 ms
// The 'good' method seems slower than the non-inverse Broyden's method, whic does not make senseto me, so I assume there is a bug
Eigen::VectorXf broyden_inv_solve(
    std::vector<std::vector<JART_VCM_v1b_var>> RRAM,
    Eigen::VectorXf Vguess, Eigen::SparseMatrix<float> G_ABCD,
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
    const float Rswl1, const float Rswl2, const float Rsbl1, const float Rsbl2,
    const float Rwl, const float Rbl,
    const bool print = false
) {
    if (RRAM.size() == 0) { return Eigen::VectorXf(0); }
    int M = RRAM.size();
    int N = RRAM[0].size();

    // Determine initial G
    Eigen::MatrixXf G(M, N);
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            float v = Vguess(i*N + j) - Vguess(i*N + j + M*N);
            G(i, j) = 1./RRAM[i][j].getResistance(v);
        }
    }

    // Calculate initial Vout
    Eigen::VectorXf Vout = solve_cam(G, Vguess, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);
    Eigen::VectorXf Fv = Vout - Vguess;

    // Calculate initial inverse Jacobian (Scaled identity matrix, finite difference, or other)
    Eigen::MatrixXf B = Eigen::MatrixXf::Identity(2*M*N, 2*M*N);

    float a = 1.;

    int it_max = 100;
    int it = 0;
    while (true) {
        if (print) { std::cout << "Norm: " << Fv.norm() << std::endl; }
        if (std::isnan(Fv.norm())) { assert(false); }
        // Check for convergence
        if (Fv.norm() < 1e-6 || it >= it_max) {
            if (print) {
                std::cout << "V:\n" << Vout << std::endl << std::endl;
                std::cout << "G:\n" << G << std::endl << std::endl;
                std::cout << "R:\n" << G.array().inverse() << std::endl << std::endl;
                std::cout << "Norm: " << (Vout - Vguess).norm() << std::endl;
                if (it >= it_max) {
                    std::cout << "Iteration limit reached: " << it << std::endl;
                } else {
                    std::cout << "Solved in " << it << " iterations" << std::endl;
                }
            }
            norm += Fv.norm();
            return Vout;
        }

        // Calculate dV
        Eigen::VectorXf dV = -B * Fv;

        // Updage Vguess
        Vguess += a * dV;
        
        // Determine G
        Eigen::MatrixXf G(M, N);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                float v = Vguess(i*N + j) - Vguess(i*N + j + M*N);
                G(i, j) = 1./RRAM[i][j].getResistance(v);
            }
        }

        // Calculate V
        Vout = solve_cam(G, Vguess, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);
        Eigen::VectorXf Fv_new = Vout - Vguess;
        Eigen::VectorXf dF = Fv_new - Fv;

        // Update inverse Jacobian
        // B += ((dV - B * dF) / (dV.transpose() * B * dF + 1e-12)) * dV.transpose() * B;
        B += ((dV - B * dF) / (dF.squaredNorm() + 1e-12)) * dF.transpose();

        Fv = Fv_new;

        it += 1;
    }
}

// Often not able to reach criterion
// Measurements are from 1 run without optimization
// For 3x3: 1.05438e-05 norm, 47 ms
// For 8x8: 4.41772e-05 norm, 766 ms
// For 16x16: 0.000429604 norm, 22821 ms
// For 32x32: 0.000861001 norm, 1295795 ms
Eigen::VectorXf broyden_solve(
    std::vector<std::vector<JART_VCM_v1b_var>> RRAM,
    Eigen::VectorXf Vguess, Eigen::SparseMatrix<float> G_ABCD,
    const Eigen::VectorXf& Vappwl1, const Eigen::VectorXf& Vappwl2,
    const Eigen::VectorXf& Vappbl1, const Eigen::VectorXf& Vappbl2,
    const float Rswl1, const float Rswl2, const float Rsbl1, const float Rsbl2,
    const float Rwl, const float Rbl,
    const bool print = false
) {
    if (RRAM.size() == 0) { return Eigen::VectorXf(0); }
    int M = RRAM.size();
    int N = RRAM[0].size();

    // Determine initial G
    Eigen::MatrixXf G(M, N);
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            float v = Vguess(i*N + j) - Vguess(i*N + j + M*N);
            G(i, j) = 1./RRAM[i][j].getResistance(v);
        }
    }

    // Calculate initial Vout
    Eigen::VectorXf Vout = solve_cam(G, Vguess, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);
    Eigen::VectorXf Fv = Vout - Vguess;

    // Calculate initial Jacobian (Scaled identity matrix, finite difference, or other)
    Eigen::MatrixXf J = Eigen::MatrixXf::Identity(2*M*N, 2*M*N);

    float a = 1.;

    int it_max = 100;
    int it = 0;
    while (true) {
        if (print) { std::cout << "Norm: " << Fv.norm() << std::endl; }
        if (std::isnan(Fv.norm())) { assert(false); }
        // Check for convergence
        if (Fv.norm() < 1e-6 || it >= it_max) {
            if (print) {
                std::cout << "V:\n" << Vout << std::endl << std::endl;
                std::cout << "G:\n" << G << std::endl << std::endl;
                std::cout << "R:\n" << G.array().inverse() << std::endl << std::endl;
                std::cout << "Norm: " << (Vout - Vguess).norm() << std::endl;
                if (it >= it_max) {
                    std::cout << "Iteration limit reached: " << it << std::endl;
                } else {
                    std::cout << "Solved in " << it << " iterations" << std::endl;
                }
            }
            norm += Fv.norm();
            return Vout;
        }

        // Calculate dV
        Eigen::VectorXf dV = J.partialPivLu().solve(-Fv);

        // Updage Vguess
        Vguess += a * dV;
        
        // Determine G
        Eigen::MatrixXf G(M, N);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                float v = Vguess(i*N + j) - Vguess(i*N + j + M*N);
                G(i, j) = 1./RRAM[i][j].getResistance(v);
            }
        }

        // Calculate V
        Vout = solve_cam(G, Vguess, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);
        Eigen::VectorXf Fv_new = Vout - Vguess;

        // Update Jacobian
        J += (((Fv_new - Fv) - J * dV) / (dV.squaredNorm() + 1e-12)) * dV.transpose();

        Fv = Fv_new;

        it += 1;
    }
}

// Often not able to reach criterion
// Measurements are from 1 run without optimization
// For 3x3: 1.95541e-05 norm, 191 ms
// For 8x8: 0.000301343 norm, 19377 ms
// For 16x16: 0.00281188 norm, 380014 ms (Seems very unstable)
Eigen::VectorXf newton_raphson_solve(
    std::vector<std::vector<JART_VCM_v1b_var>> RRAM,
    Eigen::VectorXf Vguess, Eigen::SparseMatrix<float> G_ABCD,
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

    float a = 1.;

    int it_max = 100;
    int it = 0;
    while(true) {
        // Determine G
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                float v = Vguess(i*N + j) - Vguess(i*N + j + M*N);
                G(i, j) = 1./RRAM[i][j].getResistance(v);
            }
        }

        // Calculate Vout
        Eigen::VectorXf Vout = solve_cam(G, Vguess, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);
        Eigen::VectorXf Fv = Vout - Vguess;

        if (print) { std::cout << "Norm: " << Fv.norm() << std::endl; }
        if (std::isnan(Fv.norm())) { assert(false); }
        // Check convergence
        if (Fv.norm() < 1e-6 || it >= it_max) {
            if (print) {
                std::cout << "V:\n" << Vout << std::endl << std::endl;
                std::cout << "G:\n" << G << std::endl << std::endl;
                std::cout << "R:\n" << G.array().inverse() << std::endl << std::endl;
                std::cout << "Norm: " << (Vout - Vguess).norm() << std::endl;
                if (it >= it_max) {
                    std::cout << "Iteration limit reached: " << it << std::endl;
                } else {
                    std::cout << "Solved in " << it << " iterations" << std::endl;
                }
            }
            norm += Fv.norm();
            return Vout;
        }

        // Calculate Jacobian
        Eigen::MatrixXf J = Eigen::MatrixXf::Zero(2*M*N, 2*M*N);
        for (int col = 0; col < 2*M*N; col++) {
            float delta = 1e-3;
            Vguess(col) += delta;

            for (int i = 0; i < M; i++) {
                for (int j = 0; j < N; j++) {
                    double v = Vguess(i*N + j) - Vguess(i*N + j + M*N);
                    G(i, j) = 1./RRAM[i][j].getResistance(v);
                }
            }

            Eigen::VectorXf Fj = solve_cam(G, Vguess, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl) - Vguess;

            J.col(col) = (Fj - Fv) / delta;

            Vguess(col) -= delta;
        }

        // Calculate dV
        Eigen::VectorXf dV = J.partialPivLu().solve(-Fv);

        // Update Vguess
        Vguess += a * dV;

        it++;
    }
}

// Tested up to 128 x 128
// Often not able to reach criterion
// Measurements are from 1 run without optimization
// For 3x3: solved in 10 it, 7 ms
// For 8x8: 0.000178692 norm, 505 ms (Sometimes solved in 10-12 iterations)
// For 16x16: 0.000890559 norm, 1356 ms (Very rarely returns NAN with O3 optimization)
// For 32x32: 0.00446302 norm, 4439 ms (Very rarely returns NAN with O3 optimization)
// For 64x64: 0.0255631 norm, 24641 ms
// For 128x128: 0.00284336 norm, 136759 ms (Often stalls or returns NANs)
Eigen::VectorXf fixedpoint_solve(
    std::vector<std::vector<JART_VCM_v1b_var>> RRAM,
    Eigen::VectorXf Vguess, Eigen::SparseMatrix<float> G_ABCD,
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

    float a = 0.5;

    int it_max = 100;
    int it = 0;
    while (true) {
        // Determine G
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                float v = Vguess(i*N + j) - Vguess(i*N + j + M*N);
                // if (print) {
                //     std::cout << i << " " << j << " " << v << std::endl;
                // }
                G(i, j) = (float) 1./RRAM[i][j].getResistance(v);
            }
        }
        
        // Calculate Vout
        Eigen::VectorXf Vout = solve_cam(G, Vguess, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl);
        Eigen::VectorXf Fv = Vout - Vguess;

        if (std::isnan(Fv.norm())) {
            std::cout << Vguess << std::endl;
            std::cout << G << std::endl;
            std::cout << Vappwl1 << std::endl;
            
            // std::ofstream outFile("out.txt");

            // if (!outFile) {
            //     std::cout << "No out file" << std::endl;
            //     assert(false);
            // }

            // outFile << M << std::endl;
            // outFile << N << std::endl << std::endl;
            // outFile << G << std::endl << std::endl;
            // outFile << Vguess << std::endl << std::endl;
            // outFile << Vappwl1 << std::endl << std::endl;
            // outFile << Vappwl2 << std::endl << std::endl;
            // outFile << Vappbl1 << std::endl << std::endl;
            // outFile << Vappbl2 << std::endl << std::endl;
            // outFile << Rswl1 << std::endl;
            // outFile << Rswl2 << std::endl;
            // outFile << Rsbl1 << std::endl;
            // outFile << Rsbl2 << std::endl;
            // outFile << Rwl << std::endl;
            // outFile << Rbl << std::endl;

            // outFile.close();

            assert(false);
        }

        if (print) { std::cout << "Norm: " << Fv.norm() << std::endl; }
        // Check convergence
        if (Fv.norm() < 1e-6 | it >= it_max) {
            if (print) {
                std::cout << "V:\n" << Vout << std::endl << std::endl;
                std::cout << "G:\n" << G << std::endl << std::endl;
                std::cout << "R:\n" << G.array().inverse() << std::endl << std::endl;
                std::cout << "Norm: " << Fv.norm() << std::endl;
                if (it >= it_max) {
                    std::cout << "Iteration limit reached: " << it << std::endl;
                } else {
                    std::cout << "Solved in " << it << " iterations" << std::endl;
                }
            }
            norm += Fv.norm();
            return Vout;
        }

        // Update Vguess
        Vguess += a * (Vout - Vguess);

        it++;
    }
}

int main(int argc, char* argv[]) {
    srand((unsigned int) time(0));

    int M = (argc >= 2) ? std::atoi(argv[1]) : 16;
    int N = (argc >= 3) ? std::atoi(argv[2]) : 16;

    int runs = (argc >= 4) ? std::atoi(argv[3]) : 1;

    bool print = (argc >= 5 && std::atoi(argv[4]) == 1) ? true : false;

    std::cout << "Solving " << M << "x" << N << " system in " << runs << " run(s)" << std::endl;

    long long total_time = 0;

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

    for (int i = 0; i < runs; i++) {
        float Vdd = 1.5;
        Eigen::VectorXf Vappwl1 = Eigen::VectorXf::Random(M);
        Vappwl1 = (Vappwl1.array() > 0.5).select(Eigen::VectorXf::Constant(M, Vdd), Eigen::VectorXf::Zero(M));

        // Vappwl1(0) = 0.5;
        // Vappwl1(1) = 1.;
        // Vappwl1(2) = 1.5;
        // Vappwl1 = Eigen::VectorXf::Zero(M);
        // Vappwl1(12) = Vdd;
        // Vappwl1(14) = Vdd;

        Eigen::VectorXf Vappwl2 = Eigen::VectorXf::Zero(M);
        Eigen::VectorXf Vappbl1 = Eigen::VectorXf::Zero(M);
        Eigen::VectorXf Vappbl2 = Eigen::VectorXf::Zero(M);

        if (print) {
            std::cout << "Vappwl1:\n" << Vappwl1 << std::endl << std::endl;
        }

        Eigen::VectorXf V = Eigen::VectorXf::Zero(2*M*N);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                V(i*N + j) = Vappwl1(i);
            }
        }
        auto start_time = std::chrono::high_resolution_clock::now();

        Eigen::VectorXf Vout = fixedpoint_solve(RRAM, V, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, print);
        // Eigen::VectorXf Vout = newton_raphson_solve(RRAM, V, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, print);
        // Eigen::VectorXf Vout = broyden_solve(RRAM, V, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, print);
        // Eigen::VectorXf Vout = broyden_inv_solve(RRAM, V, G_ABCD, Vappwl1, Vappwl2, Vappbl1, Vappbl2, Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl, print);

        auto end_time = std::chrono::high_resolution_clock::now();

        auto execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        std::cout << "Execution time: " << execution_time << " (ms)" << std::endl;
        total_time += execution_time;

        if (print) {
            Eigen::MatrixXf G(M, N);
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < N; j++) {
                    float v = Vout(i*N + j) - Vout(i*N + j + M*N);
                    G(i, j) = 1./RRAM[i][j].getResistance(v);
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
    
    if (runs > 1) {
        std::cout << "Average execution time: " << total_time/runs << " ms" << std::endl;
        std::cout << "Average norm: " << norm/runs << std::endl;
    }
}
