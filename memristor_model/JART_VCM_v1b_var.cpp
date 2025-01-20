#include "JART_VCM_v1b_var.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <cfloat>

// Based on this paper:
// https://ieeexplore-ieee-org.tudelft.idm.oclc.org/document/9181475

// After studying the verilog-a code the model seems to work as such:
//   The model is split into two sections which are in series:
//     The schottky tunnel section
//     The disc-plug-serial section, which can be further split into three sections, which are also in series:
//       The disc section
//       The plug section
//       The serial section
//   The voltage accross the schottky section is related to the current based on the equation given
//   Since all components are in series, the current through the schottky section is equal to the current in the discplugserial section
//   The state variable determines the resistance of the disc section
//   The resistance of the plug section is constant
//   The resistance of the serial section depends on the discplugserial current
//   The voltage across the discplugserial current then depends on the resistances and the discplugserial section
//   The discplugserial voltage and the schottky voltage determine the ion current, which in the verilog-a code is realized as a voltage source
//   The ion current determines the change in the state variable
//   Although it is not explicitly mentioned, the schottky voltage should be equal to the applied voltage minus the discplugserial voltage, since all components are in series
// The issue with this model is as follows:
//   The relationships form a circular dependency:
//     The schottky current depends on the schottky voltage
//     The discplugserial current is equal to the schottky current
//     The discplugserial voltage depends on the discplugserial current
//     The schottky voltage depends on the discplugserial voltage and the applied voltage
//   This forms a system of equations, which should be solvable by applying an initial guess for the schottky voltage and iteratively adjusting it
//   This however sounds rather slow
// Looking forward to integrating this in the crossbar model:
//   This model could be expressed as a nonlinear resistor (a resistor where the resistance depends on the voltage across it)
//   This complicates things however is still possible
//   It means the crossbar solver needs an additional iterative solver, which again sounds rather slow
//   For this to work, the model would need the following:
//     A method which applies a voltage across the device, and the model updates its state variable
//     A method which gets the resistance of the device based on an instantaneous voltage, thus does not change the state variable
//       Alternatively, a method which applies an instantaneous voltage and returns a current, but once again does not update the state variable
//     Actually, the second method is the first method with dt = 0

// Early reset bug:
//   I'm still working out why this happens. It seems to be caused by a sharp increase in I_ion, which is caused by Treal, which is caused by V_schottky/Ischottky, which might be caused by phibn
//   The point where it changes seems to also be when (V_schottky < phibn0 - phin) is triggered to be false
//   I'm yet unsure if this causes the bug or is a result of the bug
//   However, it should be noted that phibn for V_schottky=0.0799 is 0.0621074, while phibn for V_schottky = 0.08 (which triggers it condition) is 0.18, which is phibn0
//   This makes me believe that this discontinuity is the culprit of the early state transition, however I'm unsure if the discontinuity is intentional (I assume it is, since it models a diode)
//
//   After further investigation the bug arrises from the solver. Due to the discontinuity, around that point there are two possible solutions to the device state.
//   If the solver jumps over this point while searching for the solution, it might find the wrong solution, thus causing jitter and a possible early reset
//   Under real world circomstances, voltages are assumed to change continuously, so the correct solution is the one closest to the previous solution
//   This discontinuity is only present for possitive voltages and at one point, so one possible solution is to solve both sides independently and choose the right solution based on the previous solution
//   Fix jitter during reset
//
//   The split domain method has fixed the jitter and has improved the early RESET. Now it is only slightly early
//   My best guess as to what causes this slight early RESET is some error in the temperature calculation
//   A slight increase in temperature can cause it to reset slightly earlier
//   However, it is difficult to compare the temperature values around the RESET by eye. For proper comparison a simulation of the verilog-A model is needed
//   It can be noted however that the temperature peak is around 300 while it should be around 2000, so it stands to reason that the temperature around the RESET point is also slightly higher

// TODO:
//   Variability model
//   Add changing fitting parameters
//   Add changing simulation parameters (time step, adaptive time step, solve exit criterium, etc.)
//   Fix early reset
//   Fix temperature error (might also be causing early reset)
//   Fix early solve exit (f_low * f_high > 0)
//     For now I've added something that checks whether V_low and V_high are too close to each other
//     A possible reason the other check doesn't work is that there are multiple roots (specifically, an even amount of roots), so the solver wrongly things there are no solutions in the bounds

void JART_VCM_v1b_var::updateFilamentArea() {
    A = M_PI * pow(rvar, 2);
}

void JART_VCM_v1b_var::updateTemperature(double V_schottky, double V_discplugserial, double I_schottky) {
    // if (V_prev > 0.63 && V_prev < 0.642) {
    //     std::cout << Treal << ' ' << I_schottky << ' ' << Rtheff << ' ' << V_schottky << ' ' << V_discplugserial << ' ' << Rdisc << ' ' << Rplug << ' ' << Rseries << ' ' << Nreal << std::endl;
    // }
    Treal = I_schottky * (V_schottky + V_discplugserial * (Rdisc + Rplug) / (Rdisc + Rplug + Rseries)) * Rtheff + T0;
    // Treal = I_schottky * (V_schottky + I_schottky * (Rdisc + Rplug)) * Rtheff + T0;
}

double JART_VCM_v1b_var::computeSchottkyCurrent(double V_schottky, bool print) {
    double phibn;
    if (V_schottky < phibn0 - phin) {
        double psi = phibn0 - phin - V_schottky;
        phibn = phibn0 - sqrt(sqrt(pow(P_Q, 3) * zvo * Nreal * 1e26 * psi / (8 * pow(M_PI, 2) * (pow(epsphib_eff, 3)))));
        if (phibn < 0) { phibn = 0; }
    } else { phibn = phibn0; }
    phibn_out = phibn;
    // if (print) {
    //     std::cout << V_schottky << ' ' << computeSchottkyCurrent(V_schottky) << ' ' << phibn << ' ' << pow(P_Q, 3) * zvo * Nreal * 1e26 * (phibn0 - phin - V_schottky) / (8 * pow(M_PI, 2) * (pow(epsphib_eff, 3))) << std::endl;
    // }
    if (V_schottky < 0) { // TFE Schottky SET direction
        double W00 = (P_Q * P_H / (4 * M_PI)) * sqrt(zvo * Nreal * 1e26 / (mdiel * eps_eff));
        double W0 = W00 / tanh(W00 / (P_K * Treal));
        double epsprime = W00 / (W00 / (P_K * Treal) - tanh(W00 / (P_K * Treal)));
        return -A * ((Arichardson * Treal) / P_K) * sqrt(M_PI * W00 * P_Q * (fabs(V_schottky) + phibn / pow(cosh(W00/(P_K * Treal)), 2)))
        * exp(-P_Q * phibn / W0) * (exp(P_Q * fabs(V_schottky) / epsprime) - 1);
    } else { // Schottky TE RESET direction
        if (print) {
            std::cout << A * Arichardson * pow(Treal, 2) * exp(-phibn * P_Q / (P_K * Treal)) * (exp(P_Q / (P_K * Treal) * V_schottky) - 1) << ' ' << V_prev << ' ' << phibn << ' ' << V_schottky << ' ' << Nreal << std::endl;
        }
        return A * Arichardson * pow(Treal, 2) * exp(-phibn * P_Q / (P_K * Treal)) * (exp(P_Q / (P_K * Treal) * V_schottky) - 1);
    }
}

void JART_VCM_v1b_var::updateResistance(double I_discplugserial) {
    Rdisc = lvar * 1e-9 / (Nreal * 1e26 * zvo * P_Q * un * A);
    Rplug = ((lcell - lvar) * 1e-9 / (Nplug * 1e26 * zvo * P_Q * un * A));
    Rseries = RTiOx + R0 * (1 + R0 * alphaline * pow(I_discplugserial, 2) * Rthline);
}

void JART_VCM_v1b_var::updateConcentration(double I_ion, double dt) {
    double Nchange = (-trig / (A * lvar * 1e-9 * P_Q * zvo) * I_ion / 1e26) * dt;
    // Nreal = Ninitreal + Nchange;
    Nreal += Nchange;
    if (Nreal > Ndiscmax) { Nreal = Ndiscmax; }
    else if (Nreal < Ndiscmin) { Nreal = Ndiscmin; }
}

double JART_VCM_v1b_var::computeIonCurrent(double V_applied, double V_schottky, double V_discplugserial) {
    if ((Nreal < Ndiscmin && V_applied > 0) | (Nreal > Ndiscmax && V_applied < 0)) { // Keep concentration Nreal in the borders of Ndiscmin and Ndiscmax
        trig = 0;
        return 0;
    } else {
        double cvo = (Nplug + Nreal) / 2. * 1e26;
        double E_ion;
        double Flim;
        if (V_applied > 0) {
            E_ion = (V_schottky + V_discplugserial * (Rdisc + Rplug) / (Rdisc + Rplug + Rseries)) / (lcell * 1e-9);
            Rtheff = Rth0 * pow(rdet/rvar, 2) * Rtheff_scaling;
            Flim = 1 - pow(Ndiscmin/Nreal, 10);
        } else {
            E_ion = V_discplugserial * Rdisc / (Rdisc + Rplug + Rseries) / (lvar * 1e-9);
            Rtheff = Rth0 * pow(rdet/rvar, 2);
            Flim = 1 - pow(Nreal/Ndiscmax, 10);
        }
        double gamma = zvo * P_Q * E_ion * a / (M_PI * dWa * P_Q);
        double dWamin = dWa * P_Q * (sqrt(1 - pow(gamma, 2)) - gamma * M_PI/2 + gamma * asin(gamma));
        double dWamax = dWa * P_Q * (sqrt(1 - pow(gamma, 2)) + gamma * M_PI/2 + gamma * asin(gamma));
        // if (V_applied > 0.62 && V_applied < 0.72) {
        //     std::cout << V_applied << ' ' << Treal << ' ' << dWamin << ' ' << dWamax << ' ' << ' ' << zvo * P_Q * cvo * a * nyo * A * (exp(-dWamin / (P_K * Treal)) - exp(-dWamax / (P_K * Treal))) * Flim << ' ' << Nreal << std::endl;
        // }
        return zvo * P_Q * cvo * a * nyo * A * (exp(-dWamin / (P_K * Treal)) - exp(-dWamax / (P_K * Treal))) * Flim;
    }
}

std::array<double, 3> JART_VCM_v1b_var::solve_system(double V_low, double V_high, double V_applied) {
    assert(!std::isnan(V_low));
    assert(!std::isnan(V_high));
    assert(!std::isnan(V_applied));
    assert(!std::isinf(V_low));
    assert(!std::isinf(V_high));
    assert(!std::isinf(V_applied));
    if (V_low > V_high) { return {NAN, NAN, NAN}; }

    double V_schottky;
    double I_schottky;
    double V_discplugserial;

    // I_schottky = computeSchottkyCurrent(V_low);
    // updateResistance(I_schottky);
    // V_discplugserial = (Rdisc + Rplug + Rseries) * I_schottky;
    // double f_low = V_applied - V_low - V_discplugserial;
    
    // I_schottky = computeSchottkyCurrent(V_high);
    // updateResistance(I_schottky);
    // V_discplugserial = (Rdisc + Rplug + Rseries) * I_schottky;
    // double f_high = V_applied - V_high - V_discplugserial;
    // if (f_low * f_high > 0) {
    //         // std::cout << "No solution in bounds" << std::endl;
    //         // std::cout << V_applied << std::endl;
    //         return {NAN, NAN, NAN};
    // }
    
    int i = 0;
    double a = 0.5;
    while (true) {
        V_schottky = V_low * a + V_high * (1. - a);

        I_schottky = computeSchottkyCurrent(V_schottky);
        updateResistance(I_schottky);
        V_discplugserial = (Rdisc + Rplug + Rseries) * I_schottky;
        if (std::isinf(V_discplugserial)) { V_discplugserial = FLT_MAX; }
        double err = V_applied - V_discplugserial - V_schottky;
        if (fabs(err) < 1e-6) { return {V_schottky, V_discplugserial, I_schottky}; }
        // if (V_applied < 0) {
        //     if (err < 0) { V_low = V_schottky; }
        //     else { V_high = V_schottky; }
        // } else {
        //     if (err > 0) { V_low = V_schottky; }
        //     else { V_high = V_schottky; }
        // }
        if (err > 0) { V_low = V_schottky; }
        else { V_high = V_schottky; }
        
        if (i > 1e3) {
            std::cout << "Iteration limit reached, err: " << err << std::endl;
            // std::cout << V_low << " " << V_schottky << " " << V_high << std::endl;
            // std::cout << I_schottky << " " << V_discplugserial << " " << V_applied - V_discplugserial - V_schottky << std::endl;
            // return {NAN, NAN, NAN};
            return {V_schottky, V_discplugserial, I_schottky};
        }
        if (abs(V_low - V_high) < 1e-9) {
            return {NAN, NAN, NAN};
        }
        if (std::isinf(V_low) || std::isinf(V_schottky) || std::isinf(V_high)) {
            std::cout << "inf detected" << std::endl;
            // std::cout << V_low << " " << V_schottky << " " << V_high << std::endl;
            // std::cout << I_schottky << " " << V_discplugserial << " " << V_applied - V_discplugserial - V_schottky << std::endl;
            assert(false);
            return {NAN, NAN, NAN};
        }
        if (std::isnan(V_low) || std::isnan(V_schottky) || std::isnan(V_high)) {
            std::cout << "NAN detected" << std::endl;
            assert(false);
            return {NAN, NAN, NAN};
        }
        if (std::isnan(I_schottky) || std::isnan(V_discplugserial) || std::isnan(err)) {
            // std::cout << "nan detected" << std::endl;
            // std::cout << V_low << " " << V_schottky << " " << V_high << std::endl;
            // std::cout << I_schottky << " " << V_discplugserial << " " << V_applied - V_discplugserial - V_schottky << std::endl;
            return {NAN, NAN, NAN};
        }
        
        // if (V_applied > tresh) {
        //     std::cout << "V applied: " << V_applied << std::endl;
        //     std::cout << "V schottky: " << V_schottky << std::endl;
        //     std::cout << "I schottky: " << I_schottky << std::endl;
        //     std::cout << "V discplugserial: " << V_discplugserial << std::endl;
        //     std::cout << "Criterion: " << fabs((V_applied - V_discplugserial) - V_schottky) << std::endl;
        //     std::cout << "V_applied - V_discplugserial: " << V_applied - V_discplugserial << std::endl;
        //     std::cout << std::endl;
        // }
        
        i += 1;
    }
}

double JART_VCM_v1b_var::apply_voltage(double V_applied, double dt) {
    // Check voltage crossings
    if ((V_prev > -1.5e-5 && V_applied <= -1.5e-5) || (V_prev < 1.5e-5 && V_applied >= 1.5e-5)) {
        rold = rvar;
        lold = lvar;
        Nold = Nreal;
        trig = 1;
    }
    // double V_low = 0;
    // double V_high = V_applied;
    double V_schottky;
    double V_discplugserial;
    double I_schottky;
    if (abs(V_applied - V_prev) == 1e-6) {
        V_schottky = solve_prev[0];
        V_discplugserial = solve_prev[1];
        I_schottky = solve_prev[2];
    } else {
        if (V_applied < 0) {
            auto result = solve_system(V_applied, 0, V_applied);
            V_schottky = result[0];
            V_discplugserial = result[1];
            I_schottky = result[2];
        } else if (V_applied > 0) {
            auto result1 = solve_system(0, phibn0 - phin, V_applied);
            auto result2 = solve_system(phibn0 - phin, V_applied, V_applied);

            if (std::isnan(result1[0]) || std::isnan(result1[1]) || std::isnan(result1[2]) || std::isinf(result1[0]) || std::isinf(result1[1]) || std::isinf(result1[2])) {
                V_solve_bottom = 0;
            } else { V_solve_bottom = result1[0]; }
            if (std::isnan(result2[0]) || std::isnan(result2[1]) || std::isnan(result2[2]) || std::isinf(result2[0]) || std::isinf(result2[1]) || std::isinf(result2[2])) {
                V_solve_top = 0;
            } else { V_solve_top = result2[0]; }

            if (V_applied < V_prev) {
                std::swap(result1, result2);
            }
            if (std::isnan(result1[0]) || std::isnan(result1[1]) || std::isnan(result1[2]) || std::isinf(result1[0]) || std::isinf(result1[1]) || std::isinf(result1[2])) {
                V_schottky = result2[0];
                V_discplugserial = result2[1];
                I_schottky = result2[2];
            } else {
                V_schottky = result1[0];
                V_discplugserial = result1[1];
                I_schottky = result1[2];
            }
            // auto result = solve_system(0, V_applied, V_applied);
            // V_schottky = result[0];
            // V_discplugserial = result[1];
            // I_schottky = result[2];
        } else {
            V_schottky = 0;
            V_discplugserial = 0;
            I_schottky = 0;
        }
    }
    if (std::isinf(V_schottky) || std::isinf(I_schottky) || std::isinf(V_discplugserial)) {
        std::cout << "inf detected" << std::endl;
        // std::cout << "Vapplied: " << V_applied << std::endl;
        // std::cout << "Vschottky: " << V_schottky << std::endl;
        // std::cout << "Vdiscplugser: " << V_discplugserial << std::endl;
        // std::cout << "I: " << I_schottky << std::endl;
        assert(false);
    }
    if (std::isnan(V_schottky) || std::isnan(I_schottky) || std::isnan(V_discplugserial)) {
        std::cout << "nan detected" << std::endl;
        // std::cout << "Vapplied: " << V_applied << std::endl;
        // std::cout << "Vschottky: " << V_schottky << std::endl;
        // std::cout << "Vdiscplugser: " << V_discplugserial << std::endl;
        // std::cout << "I: " << I_schottky << std::endl;
        assert(false);
    }
    
    V_prev = V_applied;
    V_schottky_prev = V_schottky;
    solve_prev = {V_schottky, V_discplugserial, I_schottky};
    // f(V) = V_applied - V_schottky - V_discplugserial
    // For positive voltage, split into two regions based on shottky diode (threshold = phibn0 - phin)
    // Check if a solution exists in that range (if no solution exists f(v_low) * f(v_high) > 0)
    // If a solution exists for both, search the space on the side of the previous voltage
    // if (V_prev > 0.157 && V_prev < 0.162) {
    //     double tmp = computeSchottkyCurrent(V_schottky, true);
    //     // std::cout << V_applied << ' ' << V_schottky << ' ' << tmp << ' ' << Nreal << std::endl;
    // }
    if (dt == 0) { return I_schottky; }
    double I_ion = computeIonCurrent(V_applied, V_schottky, V_discplugserial);
    // if (V_applied > tresh && V_applied < tresh2) {
    //     std::cout << "dt: " << dt << ", V_applied: " << V_applied << ", I_ion: " << I_ion << ", Nreal: " << Nreal << std::endl;
    // }
    double N_before = Nreal;
    updateConcentration(I_ion, dt);
    // Force smaller time steps during abrupt switching
    if (fabs(Nreal - N_before) > 1e-1) {
        if (dt < 1e-9) { return I_schottky; }
        Nreal = N_before;
        for (int i = 0; i < 10; i++) {
            apply_voltage(V_applied, dt/10.);
        }
    } else {
        updateTemperature(V_schottky, V_discplugserial, I_schottky);
    }
    return I_schottky;
}

double JART_VCM_v1b_var::getResistance(double V_applied) {
    if (abs(V_applied) < 1e-6) { return Rdisc + Rplug + RTiOx; }
    else { return V_applied / apply_voltage(V_applied, 0); }
}
