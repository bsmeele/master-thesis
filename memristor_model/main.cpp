#include "JART_VCM_v1b_var.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <chrono>

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

    std::ofstream outfile("out.txt");

    if (!outfile) {
        std::cout << "No out file" << std::endl;
        return 1;
    }

    JART_VCM_v1b_var memristor = JART_VCM_v1b_var();

    const double dt = 1e-3;

    std::vector<std::array<double, 2>> Vwave;
    Vwave.push_back({0, 0});
    Vwave.push_back({-1.5, 1.5});
    Vwave.push_back({0, 3});
    Vwave.push_back({1.5, 4.5});
    Vwave.push_back({0, 6});
    Vwave.push_back({-1.5, 7.5});
    Vwave.push_back({0, 9});
    Vwave.push_back({1.5, 10.5});
    Vwave.push_back({0, 12});
    Vwave.push_back({-1.5, 13.5});
    Vwave.push_back({0, 15});
    Vwave.push_back({1.5, 16.5});
    Vwave.push_back({0, 18});

    double V = Vwave[0][0];
    double t = Vwave[0][1];

    outfile << "t V I Nreal Treal Vschottky Vdiscplugserial Rschottky Rdisc Rplug Rseries Rtotal rvar lvar" << std::endl;

    long long total_time = 0;

    for (int i = 1; i < Vwave.size(); i++) {
        double dv = (Vwave[i][0] - V) / ((Vwave[i][1] - t) / dt);
        while (t < Vwave[i][1]) {
            double I;

            auto start_time = std::chrono::high_resolution_clock::now();
            I = memristor.ApplyVoltage(V, dt);
            auto end_time = std::chrono::high_resolution_clock::now();

            auto execution_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
            total_time += execution_time;

            if (std::isnan(I)) { return 1; }
            outfile << t << " " << V << " " << I << " " << memristor.Nreal << " " << memristor.Treal
            << " " << (V - (memristor.Rdisc + memristor.Rplug + memristor.Rseries) * I) << " " << (memristor.Rdisc + memristor.Rplug + memristor.Rseries) * I
            << " " << (V - (memristor.Rdisc + memristor.Rplug + memristor.Rseries) * I)/I << " " << memristor.Rdisc << " " << memristor.Rplug << " " << memristor.Rseries
            << " " << V/I
            << " " << memristor.rvar << " " << memristor.lvar
            << std::endl;
            V += dv;
            t += dt;
        }
    }

    outfile.close();

    std::cout << "Simulated in " << total_time/1000. << " us" << std::endl;
}
