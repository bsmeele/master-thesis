#include "JART_VCM_v1b_var.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
// #include <ios>
// #include <iomanip>

int main() {
    // Applied voltage:
    //   Piecewise linear wave with:
    //     0 V at t = 0
    //     -1.5 V at t = 1.5
    //     0 V at t = 3
    //     1.5 V at t = 4.5
    //     0 V at t = 6
    // Simulation time of 8 seconds
    // Maximum step size of 3 ms

    std::ofstream outFile("out.txt");

    if (!outFile) {
        std::cout << "No out file" << std::endl;
        return 1;
    }

    JART_VCM_v1b_var memristor = JART_VCM_v1b_var();

    const double dt = 1e-3;
    // const double dt = 0.1;

    std::vector<std::array<double, 2>> V_wave;
    V_wave.push_back({0, 0});
    V_wave.push_back({-1.5, 1.5});
    V_wave.push_back({0, 3});
    V_wave.push_back({1.5, 4.5});
    V_wave.push_back({0, 6});

    double V = V_wave[0][0];
    double t = V_wave[0][1];

    // memristor.Nreal = 20;
    // double V_test = 0.160001;
    // double I_test;
    // I_test = memristor.computeSchottkyCurrent(V_test);
    // std::cout << "V: " << V_test << ", I: " << I_test << std::endl;
    // I_test = memristor.computeSchottkyCurrent(V_test-0.01);
    // std::cout << "V: " << V_test-0.01 << ", I: " << I_test << std::endl;

    outFile << "t V I Nreal Treal Rschottky Rdisc Rplug Rseries Vschottky Vdiscplugserial phibn Vsolvebottom Vsolvetop Rtotal Rtheff" << std::endl;

    for (int i = 1; i < V_wave.size(); i++) {
        // std::cout << i << std::endl;
        double dv = (V_wave[i][0] - V) / ((V_wave[i][1] - t) / dt);
        while (t < V_wave[i][1]) {
            double I;
            I = memristor.apply_voltage(V, dt);
            // std::cout << t << "s" << ": " << V << " V" << ", " << I << " A, " << "N: " << memristor.Nreal << std::endl;
            if (std::isnan(I)) { return 1; }
            outFile << t << " " << V << " " << I << " " << memristor.Nreal << " " << memristor.Treal
            << " " << (V - (memristor.Rdisc + memristor.Rplug + memristor.Rseries) * I)/I << " " << memristor.Rdisc << " " << memristor.Rplug << " " << memristor.Rseries
            << " " << (V - (memristor.Rdisc + memristor.Rplug + memristor.Rseries) * I) << " " << (memristor.Rdisc + memristor.Rplug + memristor.Rseries) * I
            << " " << memristor.phibn_out << " " << memristor.V_solve_bottom << " " << memristor.V_solve_top << " " << V/I
            << " " << memristor.Rtheff
            << std::endl;
            V += dv;
            t += dt;
        }
    }

    outFile.close();

    memristor.Treal = memristor.T0;

    std::ofstream outFile2("out2.txt");

    if (!outFile2) {
        std::cout << "No out file" << std::endl;
        return 1;
    }

    outFile2 << "V 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20" << std::endl;

    for (double v = -1.5; v <= 1.5; v += dt) {
        outFile2 << v << " ";
        for (int i = 0; i <= 20; i++) {
            memristor.Nreal = i;
            if (i == 0) { memristor.Nreal = memristor.Ndiscmin; }
            double Iout = memristor.apply_voltage(v, 0);
            if (std::abs(Iout) < 1e-12) { outFile2 << memristor.Rdisc + memristor.Rplug + memristor.Rseries; }
            else { outFile2 << v/Iout; }
            if (i < 20) { outFile2 << " "; }
        }
        outFile2 << std::endl;
    }

    outFile2.close();

    std::ofstream outFile3("out3.txt");

    if (!outFile3) {
        std::cout << "No out file" << std::endl;
        return 1;
    }

    outFile3 << "V Imin Imax V/Imin V/Imax" << std::endl;

    memristor.Treal = memristor.T0;
    for (double v = -.1; v <= .1; v += 1e-6) {
        memristor.Nreal = memristor.Ndiscmin;
        double Imin = memristor.computeSchottkyCurrent(v);
        memristor.Nreal = memristor.Ndiscmax;
        double Imax = memristor.computeSchottkyCurrent(v);
        outFile3 << v << " " << Imin << " " << Imax << " " << v/Imin << " " << v/Imax << std::endl;
    }

    outFile3.close();

    // memristor.Nreal = memristor.Ndiscmin;
    // for (double v = 0; v <= 20.; v += 1e-1) {
    //     // double I = memristor.computeSchottkyCurrent(v);
    //     // if (std::isnan(I)) {
    //     //     std::cout << v << " " << I << std::endl;
    //     // }
    //     double I = memristor.apply_voltage(v, 0);
    //     std::cout << v << " " << I << " ";
    //     std::cout << memristor.V_solve_bottom << " " << memristor.V_solve_top;
    //     std::cout << std::endl;
    // }

    // memristor.Treal = memristor.T0;
    // memristor.Nreal = memristor.Ndiscmin;
    // double Iout = memristor.apply_voltage(-1e-6, 0);
    // std::cout << Iout / -1e-6 << std::endl;
    // Iout = memristor.apply_voltage(1e-6, 0);
    // std::cout << Iout / 1e-6 << std::endl;

    // memristor.Nreal = memristor.Ndiscmax;
    // double I_test;
    // double V_test;

    // I_test = memristor.computeSchottkyCurrent(0.0799);
    // memristor.updateResistance(I_test);
    // V_test = (memristor.Rdisc + memristor.Rplug + memristor.Rseries) * I_test;
    // std::cout << 0.0799 << ' ' << I_test << ' ' << V_test << ' ' << V_test + 0.0799 << std::endl;

    // I_test = memristor.computeSchottkyCurrent(0.08);
    // memristor.updateResistance(I_test);
    // V_test = (memristor.Rdisc + memristor.Rplug + memristor.Rseries) * I_test;
    // std::cout << 0.08 << ' ' << I_test << ' ' << V_test << ' ' << V_test + 0.08 << std::endl;

    // I_test = memristor.computeSchottkyCurrent(0.0801);
    // memristor.updateResistance(I_test);
    // V_test = (memristor.Rdisc + memristor.Rplug + memristor.Rseries) * I_test;
    // std::cout << 0.0801 << ' ' << I_test << ' ' << V_test << ' ' << V_test + 0.0801 << std::endl;
    
    // for (double v = 0.07; v < 0.09; v += 0.0001) {
    //     // memristor.computeSchottkyCurrent(v, true);
    //     // std::cout << v << ' ' << memristor.computeSchottkyCurrent(v) << std::endl;
    // }

    // memristor.Nreal = memristor.Ndiscmax;
    // memristor.Treal = memristor.T0;
    // I_test = memristor.apply_voltage(0.16, 0);
    // memristor.updateResistance(I_test);
    // V_test = 0.16 - (memristor.Rdisc + memristor.Rplug + memristor.Rseries)*I_test;
    // std::cout << V_test << " " << (memristor.Rdisc + memristor.Rplug + memristor.Rseries)*I_test << std::endl;

    // std::cout << "V_schottky   V_discplugserial   V_applied   I_shottky" << std::endl;

    // memristor.Nreal = memristor.Ndiscmax;
    // memristor.Treal = memristor.T0;
    // V_test = 0.117506;
    // I_test = memristor.computeSchottkyCurrent(V_test);
    // memristor.updateResistance(I_test);
    // std::cout << V_test << "     " << (memristor.Rdisc + memristor.Rplug + memristor.Rseries)*I_test << "          " << (memristor.Rdisc + memristor.Rplug + memristor.Rseries)*I_test + V_test << "    " << I_test << std::endl;

    // memristor.Nreal = memristor.Ndiscmax;
    // memristor.Treal = memristor.T0;
    // V_test = 0.00663375;
    // I_test = memristor.computeSchottkyCurrent(V_test);
    // memristor.updateResistance(I_test);
    // std::cout << V_test << "   " << (memristor.Rdisc + memristor.Rplug + memristor.Rseries)*I_test << "           " << (memristor.Rdisc + memristor.Rplug + memristor.Rseries)*I_test + V_test << "    " << I_test << std::endl;

    // memristor.Nreal = memristor.Ndiscmax;
    // memristor.Treal = 440.;
    // V_test = 0.1481409;
    // I_test = memristor.computeSchottkyCurrent(V_test);
    // memristor.updateResistance(I_test);
    // std::cout << V_test << "     " << (memristor.Rdisc + memristor.Rplug + memristor.Rseries)*I_test << "            " << (memristor.Rdisc + memristor.Rplug + memristor.Rseries)*I_test + V_test << "    " << I_test << std::endl;

    // memristor.Nreal = memristor.Ndiscmax;
    // memristor.Treal = 440.;
    // V_test = 0.01619361;
    // I_test = memristor.computeSchottkyCurrent(V_test);
    // memristor.updateResistance(I_test);
    // std::cout << V_test << "    " << (memristor.Rdisc + memristor.Rplug + memristor.Rseries)*I_test << "           " << (memristor.Rdisc + memristor.Rplug + memristor.Rseries)*I_test + V_test << "    " << I_test << std::endl;

    // memristor.Nreal = memristor.Ndiscmax;
    // memristor.Treal = 600.;
    // V_test = 0.07999;
    // I_test = memristor.computeSchottkyCurrent(V_test);
    // memristor.updateResistance(I_test);
    // std::cout << V_test << "    " << (memristor.Rdisc + memristor.Rplug + memristor.Rseries)*I_test << "           " << (memristor.Rdisc + memristor.Rplug + memristor.Rseries)*I_test + V_test << "    " << I_test << std::endl;
}