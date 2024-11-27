#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <chrono>

// Based on:
// https://ieeexplore-ieee-org.tudelft.idm.oclc.org/document/9108292

Eigen::VectorXf solve_fcu(const Eigen::MatrixXf& G, const Eigen::VectorXf& Vin, const float Rcol, const float Rrow, const float Rsense, bool print = false) {
    // auto start_time = std::chrono::high_resolution_clock::now();

    int M = G.rows();
    int N = G.cols();

    // ----- Step 1: Formulate Column Linear Systems -----
    // A consists of N matrices, which are each MxM
    // std::vector<Eigen::MatrixXf> A;
    // for (int j = 0; j < N; j++) {
    //     Eigen::MatrixXf Aj = Eigen::MatrixXf::Zero(M, M);
    //     for (int row = 0; row < M; row++) {
    //         for (int i = 0; i < M; i++) {
    //             if (i == row) { Aj(row, i) = -1; }
    //             else if (i > row) { Aj(row, i) = (i - row) * G(i, j) * Rcol; }
    //         }
    //     }

    //     A.push_back(Aj);
    // }

    // if (print) {
    //     for (int j = 0; j < N; j++) {
    //         std::cout <<  "A" << j << std::endl << A[j] << std::endl << std::endl;
    //     }
    // }

    // // J is a column vector consisting of M elements
    // Eigen::VectorXf J(M);
    // for (int i = 0; i < M; i++) {
    //     J(i) = Rsense + (M - (i + 1)) * Rcol;
    // }

    // if (print) {
    //     std::cout << "J:\n" << J << std::endl << std::endl;
    // }

    // // K consists of N row vectors, each with M elements
    // std::vector<Eigen::MatrixXf> K;
    // for (int j = 0; j < N; j++) {
    //     Eigen::RowVectorXf Kj(M);
    //     for (int i = 0; i < M; i++) {
    //         Kj(i) = G(i, j);
    //     }
    //     K.push_back(Kj);
    // }

    // if (print) {
    //     for (int j = 0; j < N; j++) {
    //         std::cout <<  "K" << j << std::endl << K[j] << std::endl << std::endl;
    //     }
    // }

    // auto timestamp1 = std::chrono::high_resolution_clock::now();
    // auto execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp1 - start_time).count();
    // std::cout << "Setup time: " << execution_time << " (ms)" << std::endl;
    // timestamp1 = std::chrono::high_resolution_clock::now();

    // ----- Step 2: Merge Column Linear Systems -----
    // COLmat is the direct sum of (Aj + J * Kj), thus is a M*N x M*N matrix
    Eigen::MatrixXf COLmat = Eigen::MatrixXf::Zero(M*N, M*N);
    for (int j = 0; j < N; j++) {
        for (int row = 0; row < M; row++) {
            for (int i = 0; i < M; i++) {
                if (i == row) { COLmat(row + j*M, i + j*M) = 1 + (Rsense + (M - (row + 1)) * Rcol) * G(i, j); }
                else if (i > row) { COLmat(row + j*M, i + j*M) = -(i - row) * G(i, j) * Rcol + (Rsense + (M - (row + 1)) * Rcol) * G(i, j); }
                else { COLmat(row + j*M, i + j*M) = (Rsense + (M - (row + 1)) * Rcol) * G(i, j); }
            }
        }
        // Eigen::MatrixXf tmp = -A[j] + J * K[j];
        // COLmat.block(j*M, j*M, tmp.rows(), tmp.cols()) = tmp;
    }

    // auto timestamp2 = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp2 - timestamp1).count();
    // std::cout << "COLmat time: " << execution_time << " (ms)" << std::endl;
    // timestamp2 = std::chrono::high_resolution_clock::now();

    if (print) {
        std::cout << "COLmat:\n" << COLmat << std::endl << std::endl;
    }

    // Gmat is the direct sum of Kj, thus is a N x M*N matrix
    Eigen::MatrixXf Gmat = Eigen::MatrixXf::Zero(N, M*N);
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < M; i++) {
            Gmat(j, i + j*M) = G(i, j);
        }
    }

    // auto timestamp3 = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp3 - timestamp2).count();
    // std::cout << "Gmat time: " << execution_time << " (ms)" << std::endl;
    // timestamp3 = std::chrono::high_resolution_clock::now();

    if (print) {
        std::cout << "Gmat:\n" << Gmat << std::endl << std::endl;
    }

    // ----- Step 3: Formulate Row Linear Systems -----
    // B consists of M matrices, which are each NxN
    // std::vector<Eigen::MatrixXf> B;
    // for (int i = 0; i < M; i++) {
    //     Eigen::MatrixXf Bi(N, N);
    //     for (int row = 0; row < N; row++) {
    //         for (int j = 0; j < N; j++) {
    //             if (row <= j) { Bi(row, j) = G(i, j) * (row + 1); }
    //             else { Bi(row, j) = G(i, j) * (j + 1); }
    //         }
    //     }

    //     B.push_back(Bi*Rrow);
    // }

    // if (print) {
    //     for (int i = 0; i < M; i++) {
    //         std::cout <<  "B" << i << std::endl << B[i] << std::endl << std::endl;
    //     }
    // }
    
    // ----- Step 4: Merge Row Linear Systems -----
    // ROWmat is the direct sum of all Bi matrices, thus is a N*M x N*M matrix
    Eigen::MatrixXf ROWmat = Eigen::MatrixXf::Zero(N*M, N*M);
    for (int i = 0; i < M; i++) {
        for (int row = 0; row < N; row++) {
            for (int j = 0; j < N; j++) {
                if (row <= j) { ROWmat(row + i*N, j + i*N) = G(i, j) * (row + 1); }
                else { ROWmat(row + i*N, j + i*N) = G(i, j) * (j + 1); }
            }
        }
        // ROWmat.block(i*N, i*N, B[i].rows(), B[i].cols()) = B[i];
    }
    ROWmat *= Rrow;

    // auto timestamp4 = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp4 - timestamp3).count();
    // std::cout << "ROWmat time: " << execution_time << " (ms)" << std::endl;
    // timestamp4 = std::chrono::high_resolution_clock::now();

    if (print) {
        std::cout << "ROWmat:\n" << ROWmat << std::endl << std::endl;
    }

    // ----- Step 5: Eliminate Internal Variables -----
    // ROWmatA and CVrowINA are constructed by applying row swaps to ROWmat and CVrowIN respectively, thus also have the same shapes
    Eigen::MatrixXf ROWmatAint(N*M, N*M);
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < M; i++) {
            ROWmatAint.row(j*M + i) = ROWmat.row(i*N + j);
        }
    }
    Eigen::MatrixXf ROWmatA(N*M, N*M);
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < M; i++) {
            ROWmatA.col(j*M + i) = ROWmatAint.col(i*N + j);
        }
    }

    // auto timestamp5 = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp5 - timestamp4).count();
    // std::cout << "ROWmatA time: " << execution_time << " (ms)" << std::endl;
    // timestamp5 = std::chrono::high_resolution_clock::now();

    if (print) {
        std::cout << "ROWmatA:\n" << ROWmatA << std::endl << std::endl;
    }

    // NETmat is constructe as such: Gmat * (COLmat + ROWmatA)^-1, thus is a N x N*M matrix
    // Eigen::MatrixXf NETmat = Gmat * (COLmat + ROWmatA).inverse();

    // Alternative:
    // Rewrite NEtmat = Gmat * (COLmat + ROWmatA)^-1 as (COLmat + ROWmatA)^T * NETmat^T = Gmat^T
    // Then use a solver to solve for NETmat^T
    // Finally transpose to get NETmat    
    Eigen::SparseMatrix<float> COLROWmat = (COLmat + ROWmatA).sparseView();
    Eigen::MatrixXf Gmat_t = Gmat.transpose();
    Eigen::SparseMatrix<float> Gmat_t_sparse = Gmat_t.sparseView();

    // auto timestamp6 = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp6 - timestamp5).count();
    // std::cout << "Solver setup time: " << execution_time << " (ms)" << std::endl;
    // timestamp6 = std::chrono::high_resolution_clock::now();

    Eigen::BiCGSTAB<Eigen::SparseMatrix<float>> solver;

    solver.compute(COLROWmat.transpose());
    Eigen::SparseMatrix<float> NETmat = solver.solve(Gmat_t_sparse).transpose();

    // auto timestamp7 = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp7 - timestamp6).count();
    // std::cout << "Solve time: " << execution_time << " (ms)" << std::endl;
    // timestamp7 = std::chrono::high_resolution_clock::now();

    if (print) {
        std::cout << "NETmat:\n" << NETmat << std::endl << std::endl;
    }

    // ----- Step 6: Reduce Matrix Dimensions -----
    // NETmatC is constructed by using column sums of NETmat such that CVrowINA can be reduced to Vin, thus is a NxM matrix
    Eigen::MatrixXf NETmatC = Eigen::MatrixXf::Zero(N, M);
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            NETmatC.col(i) += NETmat.col(j*M + i);
        }
    }

    // auto timestamp8 = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp8 - timestamp7).count();
    // std::cout << "NETmatC time: " << execution_time << " (ms)" << std::endl;
    // timestamp8 = std::chrono::high_resolution_clock::now();

    if (print) {
        std::cout << "NETmatC:\n" << NETmatC << std::endl << std::endl;
    }

    Eigen::VectorXf Iout = NETmatC * Vin;

    // auto end_time = std::chrono::high_resolution_clock::now();
    // execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - timestamp8).count();
    // std::cout << "Iout time: " << execution_time << " (ms)" << std::endl;

    return Iout;
}

int main(int argc, char* argv[]) {
    long long total_time = 0;

    srand((unsigned int) time(0));

    int M = (argc >= 2) ? std::atoi(argv[1]) : 16;
    int N = (argc >= 3) ? std::atoi(argv[2]) : 16;

    int runs = (argc >= 4) ? std::atoi(argv[3]) : 1;

    bool print = false;
    if (std::atoi(argv[4]) == 1) { print = true; }

    for (int i = 0; i < runs; i++) {
        float Rmin = 100.;
        float Rmax = 1000.;
        float Vdd = 5.;

        float Rcol = 2.;
        float Rrow = 3.;
        float Rsense = 5.;

        if (print) {
            std::cout << "Rcol: " << Rcol << std::endl;
            std::cout << "Rrow: " << Rrow << std::endl;
            std::cout << "Rsense: " << Rsense << std::endl << std::endl;
        }

        // Construct G matrix with random values
        Eigen::MatrixXf R = Eigen::MatrixXf::Random(M, N);
        R = (R + Eigen::MatrixXf::Constant(M, N, 1.0)) * 0.5 * (Rmax - Rmin) + Eigen::MatrixXf::Constant(M, N, Rmin);
        R(0, 0) = 10.;
        R(0, 1) = 15.;
        R(0, 2) = 20.;
        R(1, 0) = 25.;
        R(1, 1) = 30.;
        R(1, 2) = 35.;
        R(2, 0) = 40.;
        R(2, 1) = 45.;
        R(2, 2) = 50.;

        if (print) {
            std::cout << "R:\n" << R << std::endl << std::endl;
        }

        Eigen::MatrixXf G = R.cwiseInverse();

        if (print) {
            std::cout << "G:\n" << G << std::endl << std::endl;
        }

        // Construct Vin with random values
        Eigen::VectorXf Vin = Eigen::VectorXf::Random(M);
        Vin = (Vin.array() > 0.5).select(Eigen::VectorXf::Constant(M, Vdd), Eigen::VectorXf::Zero(M));
        Vin(0) = 5;
        Vin(1) = 7;
        Vin(2) = 9;

        if (print) {
            std::cout << "Vin:\n" << Vin << std::endl << std::endl;
        }

        // Eigen::VectorXf Iout(N)

        auto start_time = std::chrono::high_resolution_clock::now();

        Eigen::VectorXf Iout = solve_fcu(G, Vin, Rcol, Rrow, Rsense, print);

        auto end_time = std::chrono::high_resolution_clock::now();

        if (print) {
            std::cout << Iout << std::endl << std::endl;
        }

        auto execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        
        std::cout << "Execution time: " << execution_time << " ms" << std::endl;

        total_time += execution_time;
    }

    if (runs > 1) {
        std::cout << "Average execution time: " << total_time/runs << " ms" << std::endl;
    }

    return 0;
}
