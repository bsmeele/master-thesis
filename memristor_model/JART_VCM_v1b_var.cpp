#include "JART_VCM_v1b_var.h"
#include "../simulation_settings.h"

#include <cmath>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <cfloat>

// Based on this paper:
// https://ieeexplore-ieee-org.tudelft.idm.oclc.org/document/9181475

// The majority of this code is a recreation of the Verilog-A code of the JART VCM v1b var model from this website:
// https://www.emrl.de/JART.html#Artikel_3
// As the original model was written in Verilog-A, some adjustment have been made and additional simulation tools have been developed

// Updates the filament area based on the filament radius
// Used for the variability model
// As the variability model is not implemented this function is unused
void JART_VCM_v1b_var::UpdateFilamentArea() {
    A = M_PI * pow(rvar, 2);
}

// Updated the temperature of the device based on the power dissipated by the device
void JART_VCM_v1b_var::UpdateTemperature(double V_schottky, double V_discplugserial, double I_schottky) {
    Treal = I_schottky * (V_schottky + V_discplugserial * (Rdisc + Rplug) / (Rdisc + Rplug + Rseries)) * Rtheff + T0;
}

// Computes the schottky current based on an applied schottky voltage
double JART_VCM_v1b_var::ComputeSchottkyCurrent(double V_schottky) {
    // The calculation of phibn includes some additions compared to as it's presented in the reference paper
    // These additions are based on the Verilog-A code of the device and are as followed:
    //   If phibn would be negative, it is instead taken as 0
    //   If phibn would be complex, it is instead taken as phibn0
    double phibn;
    if (V_schottky < phibn0 - phin) {
        double psi = phibn0 - phin - V_schottky;
        phibn = phibn0 - sqrt(sqrt(pow(P_Q, 3) * zvo * Nreal * 1e26 * psi / (8 * pow(M_PI, 2) * (pow(epsphib_eff, 3)))));
        if (phibn < 0) { phibn = 0; }
    } else { phibn = phibn0; }

    if (V_schottky < 0) { // TFE Schottky SET direction
        double W00 = (P_Q * P_H / (4 * M_PI)) * sqrt(zvo * Nreal * 1e26 / (mdiel * eps_eff));
        double W0 = W00 / tanh(W00 / (P_K * Treal));
        double epsprime = W00 / (W00 / (P_K * Treal) - tanh(W00 / (P_K * Treal)));
        return -A * ((Arichardson * Treal) / P_K) * sqrt(M_PI * W00 * P_Q * (fabs(V_schottky) + phibn / pow(cosh(W00/(P_K * Treal)), 2)))
        * exp(-P_Q * phibn / W0) * (exp(P_Q * fabs(V_schottky) / epsprime) - 1);
    } else { // Schottky TE RESET direction
        return A * Arichardson * pow(Treal, 2) * exp(-phibn * P_Q / (P_K * Treal)) * (exp(P_Q / (P_K * Treal) * V_schottky) - 1);
    }
}

// Updates the resistance of the various elements of the device
// based on the current through the device and the internal state of the device
void JART_VCM_v1b_var::UpdateResistance(double I_discplugserial) {
    Rdisc = lvar * 1e-9 / (Nreal * 1e26 * zvo * P_Q * un * A);
    Rplug = ((lcell - lvar) * 1e-9 / (Nplug * 1e26 * zvo * P_Q * un * A));
    Rseries = RTiOx + R0 * (1 + R0 * alphaline * pow(I_discplugserial, 2) * Rthline);
}

// Updates the internal state of the device
// Instead of integration, as it is presented in the reference paper,
// the state variable is updated by multiplying a rate of change with a delta time
void JART_VCM_v1b_var::UpdateConcentration(double I_ion, double dt) {
    double Nchange = (-trig / (A * lvar * 1e-9 * P_Q * zvo) * I_ion / 1e26) * dt;
    Nreal += Nchange;
    if (Nreal > Ndiscmax) { Nreal = Ndiscmax; }
    else if (Nreal < Ndiscmin) { Nreal = Ndiscmin; }
}

// Computes the ion current through the device
// The ion current is exclusively used in updating the state variable of the device
double JART_VCM_v1b_var::ComputeIonCurrent(double V_applied, double V_schottky, double V_discplugserial) {
    if ((Nreal <= Ndiscmin && V_applied > 0) | (Nreal >= Ndiscmax && V_applied < 0)) { // Keep concentration Nreal in the borders of Ndiscmin and Ndiscmax
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
        return zvo * P_Q * cvo * a * nyo * A * (exp(-dWamin / (P_K * Treal)) - exp(-dWamax / (P_K * Treal))) * Flim;
    }
}

// Attempts to find compute the schottky voltage, device current, and voltage accros the remainder of the device based on some applied voltage
// The solving method used in this function is a simple fixed point method
std::array<double, 3> JART_VCM_v1b_var::SolveFixedpoint(double V_guess, double V_applied) {
    assert(!std::isnan(V_guess));
    assert(!std::isinf(V_guess));
    assert(!std::isnan(V_applied));
    assert(!std::isinf(V_applied));

    double V_schottky;
    double I_schottky;
    double V_discplugserial;

    int it = 0;
    double a = 1.;
    double a_ref = 0.2;

    double Vn_min_1 = 0;
    double Vn_min_2 = 0;

    double Fv_prev = INFINITY;
    while (true) {
        I_schottky = ComputeSchottkyCurrent(V_guess);
        UpdateResistance(I_schottky);
        V_discplugserial = (Rdisc + Rplug + Rseries) * I_schottky;
        if (std::isinf(V_discplugserial)) { V_discplugserial = FLT_MAX; }

        double V_schottky = V_applied - V_discplugserial;

        double Fv = V_schottky - V_guess;
        if (fabs(Fv) < memristor_fixedpoint_criterion) { return {V_schottky, V_discplugserial, I_schottky}; }
        
        double a = (it > 100) ? 5e-2 : 5e-3;
        V_guess = V_guess + a * Fv;

        if (V_applied > 0) {
            if (V_guess > V_applied) { V_guess = V_applied; }
            else if (V_guess < 0) { V_guess = 0; }
        } else {
            if (V_guess < V_applied) { V_guess = V_applied; }
            else if (V_guess > 0) { V_guess = 0; }
        }

        if (std::isinf(V_guess)) {
            return {NAN, NAN, NAN};
        }
        if (std::isnan(V_guess)) {
            return {NAN, NAN, NAN};
        }
        if (std::isnan(I_schottky) || std::isnan(V_discplugserial) || std::isnan(Fv)) {
            return {NAN, NAN, NAN};
        }
        if (it > memristor_fixedpoint_it_max) {
            if (memristor_fixedpoint_warn_iteration_limit) {
                std::cout << "Memristor fixed-point iteration limit reached with error: " << Fv << std::endl;
            }
            return {V_schottky, V_discplugserial, I_schottky};
        }

        it += 1;
    }
}

// Calculates the current through the device and updated internal variables based on some applied voltage and some delta time
// Due to the nature of the formulation of the device (as can be seen in the reference paper),
// calculation of the voltage accross the subcomponents and current through the device requires some root finding method
// Supplying a dt of 0 prevents this function from updated ingernal variables and merely return the current through the device
double JART_VCM_v1b_var::ApplyVoltage(double V_applied, double dt) {
    // Check voltage crossings
    if (dt != 0) {
        if ((V_prev > -1.5e-5 && V_applied <= -1.5e-5) || (V_prev < 1.5e-5 && V_applied >= 1.5e-5)) {
            rold = rvar;
            lold = lvar;
            Nold = Nreal;
            trig = 1;
        }
        // Unused variation model code:
        // if (V_applied < -2e-5 && trig == 1) {  // SET at negative voltage
        //     rvar = rold + (rnew - rold) * ((Nreal - Nold) / (Ndiscmax - Nold));
        //     lvar = lold + (lnew - lold) * ((Nreal - Nold) / (Ndiscmax - Nold));
        //     UpdateFilamentArea();
        // } else if (V_applied > 2e-5 && trig == 1) {  // RESET at postivive voltage
        //     rvar = rold + (rnew - rold) * ((Nold - Nreal) / (Nold - Ndiscmin));
        //     lvar = lold + (lnew - lold) * ((Nold - Nreal) / (Nold - Ndiscmin));
        //     UpdateFilamentArea();
        // }
    }

    double V_schottky;
    double V_discplugserial;
    double I_schottky;
    std::array<double, 3> result;

    switch (memristor_solving_method) {
        case 0:
            if (V_applied < 0) {
                result = SolveBisection(V_applied, 0, V_applied);
                V_schottky = result[0];
                V_discplugserial = result[1];
                I_schottky = result[2];
            } else if (V_applied < phibn0 - phin) {
                result = SolveBisection(0, V_applied, V_applied);
                V_schottky = result[0];
                V_discplugserial = result[1];
                I_schottky = result[2];
            } else if (V_applied >= phibn0 - phin) {
                auto result1 = SolveBisection(0, phibn0 - phin, V_applied);
                auto result2 = SolveBisection(phibn0 - phin, V_applied, V_applied);
            
                if (fabs(result2[0] - V_schottky_prev) < fabs(result1[0] - V_schottky_prev)) {
                    V_schottky = result2[0];
                    V_discplugserial = result2[1];
                    I_schottky = result2[2];
                } else {
                    V_schottky = result1[0];
                    V_discplugserial = result1[1];
                    I_schottky = result1[2];
                }
            } else {
                V_schottky = 0;
                V_discplugserial = 0;
                I_schottky = 0;
            }
            break;
        case 1:
            result = SolveFixedpoint(V_schottky_prev, V_applied);
            V_schottky = result[0];
            V_discplugserial = result[1];
            I_schottky = result[2];
            break;
        default:
            result = SolveFixedpoint(V_schottky_prev, V_applied);
            V_schottky = result[0];
            V_discplugserial = result[1];
            I_schottky = result[2];
    }

    if (std::isinf(V_schottky) || std::isinf(I_schottky) || std::isinf(V_discplugserial)) {
        return NAN;
    }
    if (std::isnan(V_schottky) || std::isnan(I_schottky) || std::isnan(V_discplugserial)) {
        return NAN;
    }
    
    if (dt != 0) {
        V_prev = V_applied;
        V_schottky_prev = V_schottky;
    } else {
        return I_schottky;
    }

    double I_ion = ComputeIonCurrent(V_applied, V_schottky, V_discplugserial);
    double N_before = Nreal;
    UpdateConcentration(I_ion, dt);

    // Force smaller time steps during abrupt switching
    if (fabs(Nreal - N_before) > memristor_dynamic_time_step_N_limit) {
        if (dt < memristor_dynamic_time_step_t_limit) { return I_schottky; }
        Nreal = N_before;
        for (int i = 0; i < memristor_dynamic_time_step_time_division; i++) {
            ApplyVoltage(V_applied, dt/memristor_dynamic_time_step_time_division);
        }
    } else {
        UpdateTemperature(V_schottky, V_discplugserial, I_schottky);
    }
    return I_schottky;
}

// Computes the resistance of the device by calculating the current for some applied voltage
// To avoid instability, this function returns some default value (based on internal variables) for low voltages
double JART_VCM_v1b_var::GetResistance(double V_applied) {
    if (fabs(V_applied) < memristor_get_resistance_voltage_threshold) {
        return Rdisc + Rplug + RTiOx + R0;
    }
    else { return V_applied / ApplyVoltage(V_applied, 0); }
}

void JART_VCM_v1b_var::SetWeight(bool weight) {
    if (weight) { Nreal = Ndiscmax; }
    else { Nreal = Ndiscmin; }
    // else { Nreal = Ninit; }
    Treal = T0;
}
