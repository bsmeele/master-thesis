#include "JART_VCM_v1b_var.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>

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

    std::vector<std::array<double, 2>> V_wave;
    V_wave.push_back({0, 0});
    V_wave.push_back({-1.5, 1.5});
    V_wave.push_back({0, 3});
    V_wave.push_back({1.5, 4.5});
    V_wave.push_back({0, 6});

    double V = V_wave[0][0];
    double t = V_wave[0][1];

    outFile << "t V I Nreal Treal Vschottky Vdiscplugserial Rschottky Rdisc Rplug Rseries Rtotal" << std::endl;

    for (int i = 1; i < V_wave.size(); i++) {
        double dv = (V_wave[i][0] - V) / ((V_wave[i][1] - t) / dt);
        while (t < V_wave[i][1]) {
            double I;
            I = memristor.apply_voltage(V, dt);
            if (std::isnan(I)) { return 1; }
            outFile << t << " " << V << " " << I << " " << memristor.Nreal << " " << memristor.Treal
            << " " << (V - (memristor.Rdisc + memristor.Rplug + memristor.Rseries) * I) << " " << (memristor.Rdisc + memristor.Rplug + memristor.Rseries) * I
            << " " << (V - (memristor.Rdisc + memristor.Rplug + memristor.Rseries) * I)/I << " " << memristor.Rdisc << " " << memristor.Rplug << " " << memristor.Rseries
            << " " << V/I
            << std::endl;
            V += dv;
            t += dt;
        }
    }

    outFile.close();
}
