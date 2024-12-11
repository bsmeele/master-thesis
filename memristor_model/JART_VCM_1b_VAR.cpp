#include <iostream>
#include <cmath>
#include <algorithm>

// Constants
const float M_PI = 3.1415927;    // Pi
const float  P_Q = 1.6022e-19;   // Elementary charge [C]
const float P_K = 1.38065e-23;   // Boltzman constant [J/K]
const float P_EPS0 = 8.6549e-12; // Permittivity of a vacuum [F/m]

const float P_H = 6.626e-34;     // Planck constant [Js]

class JART_VCM_1b_VAR {

    // Pyisical constants do not change!
    const float Arichardson = 6.01e5; // Richardson's constant [A/m^2K^2]
    const float mdiel = 9.10938e-31;  // electron rest mass [kg]
    const float zvo = 2;              // oxygen vacancy charge number
    const float T0 = 293;             // ambient temperature [K]

    // Fitting parameters
    float eps = 17;              // from [10:25], static hafnium oxide permittivity
    float epsphib = 5.5;         // hafnium oxide permittivity related to image force barrier lowering
    float phibn0 = 0.18;         // from [0.1:0.5], nominal schottky barrier height [eV]
    float phin = 0.1;            // from [0.1:0.3], energy level difference between the Fermi level in the oxide and the oxide conduction band edge [eV]
    float un = 4e-6;             // from [1e-6:1e-5], electron mobility [m^2/Vs]
    float Ndiscmax = 20;         // from [0.001:1100], maximum oxygen vacancy concentration in the disc [10^26/m^3]
    float Ndiscmin = 0.008;      // from [0.0001:100], minimum oxygen vacancy concentration in the disc [10^26/m^3]
    float Ninit = 0.008;         // from [0.0001:1000], initial oxygen vacancy concentration in the disc [10^26/m^3]
    float Nplug = 20;            // from [0.001:100], oxygen vacancy concentration in the plug [10^26/m^3]
    float a = 0.25e-9;           // from [0.1e-9:1e-9], ion hopping distance [m]
    float nyo = 2e13;            // from [1e10:1e14], attempt frequency [Hz]
    float dWa = 1.35;            // from [0.8:1.5], activation energy [eV]
    float Rth0 = 15.72e6;        // from [1e6:20e6], thermal resistance of the Hafnium Oxide [K/W]
    float rdet = 45e-9;          // from [5e-9:100e-9], radius of the filament [m]
    float rnew = 45e-9;          // from [5e-9:100e-9], radius of the filament [m]
    float lcell = 3;             // from [2:5], length of disc and plug region [nm]
    float ldet = 0.4;            // from [0.1:5], length of the disc region [nm]
    float lnew = 0.4;            // from [0.1:5], length of the disc region [nm]
    float Rtheff_scaling = 0.27; // from [0.1:1], scaling factor for RESET
    float RTiOx = 650;           // from [0:5000], series resistance of the TiOx layer [Ohm]
    float R0 = 719.2437;         // Resistance at T0 [Ohm]
    float Rthline = 90471.47;    // thermal conductivity of the Platinum and Titanium [W/mK]
    float alphaline = 3.92e-3;   // temperature coefficient [1/K]

    // Variables
    float A;
    float Flim;
    float phibn;
    float Rtheff;
    float psi;
    float W00;
    float W0;
    float epsprime;
    float Ischottkytunnel;
    float cvo;
    float Treal;
    float Nreal;
    float Nchange;
    float E_ion;
    float dWamin;
    float dWamax;
    float gamma;
    float Rdisc;
    float Rplug;
    float Rseries;
    float rvar;
    float lvar;
    float lold;
    float Nold;
    float rold;

    float Ninitreal;
    int trig;

    float eps_eff;     // static hafnium oxide permittivity
    float epsphib_eff; // hafnium oxide permittivity related to image force barrier lowering

    // Initial step
    JART_VCM_1b_VAR() : Ninitreal(Ninit), trig(1), rvar(rnew), lvar(lnew), lold(lnew), rold(rnew), Nold(Ninit) {
        eps_eff = eps * P_EPS0;
        epsphib_eff = epsphib * P_EPS0;
    }

    // Variability model skipped for simplicity, to be added later

    void updateFilamentArea() {
        A = M_PI * pow(rvar, 2);
    }

    // Something about forcing a small time step during abrupt switching
    // if ((abs(V(slopeN,gnd))>1e-8))
	// begin
	// 	$bound_step(5e-12);
	// end

    void computeTemperature(float I_discplugseries, float V_disc, float V_plug, float V_schottky) {
        Treal = T0 + I_discplugseries * (V_disc + V_plug + V_schottky) * Rtheff;
    }

    float computeSchottkyCurrent(float V_Schottky) {
        if (V_Schottky < phibn0 - phin) {
            psi = phibn0 - phin - V_Schottky;
            phibn = phibn0 - sqrt(sqrt(pow(P_Q, 3) * zvo * Nreal * 1e26 * psi / (8 * pow(M_PI,2) * (pow(epsphib_eff,3)))));
            if (phibn < 0) { phibn = 0; }
        } else {
            psi = 0;
            phibn = phibn0;
        }

        if (V_Schottky < 0) { // TFE Schottky SET direction
            W00 = (P_Q * P_H / (4 * M_PI)) * sqrt(zvo * Nreal * 1e26 / (mdiel * eps_eff));
            W0 = W00 / tanh(W00 / (P_K * Treal));
            epsprime = W00 / (W00 / (P_K * Treal) - tanh(W00 / (P_K * Treal)));
            return -A * ((Arichardson * Treal) / P_K) * sqrt(M_PI * W00 * P_Q * (abs(V_Schottky) + phibn / pow(cosh(W00/(P_K * Treal)),2)))
                * exp(-P_Q * phibn / W0) * (exp(P_Q * abs(V_Schottky) / epsprime) - 1);
        } else { // Schotky TE RESET direction
            return A * Arichardson * pow(Treal, 2) * exp(-phibn * P_Q / (P_K * Treal)) * (exp(P_Q / (P_K * Treal) * V_Schottky) - 1);
        }
    }

    float computeDiscplugserialVoltage(float I_discplugserial) {
        Rdisc = lvar * 1e-9 / (Nreal * 1e26 * zvo * P_Q * un * A);
        Rplug = ((lcell - lvar) * 1e-9 * (Nplug * 1e26 * zvo * P_Q * un * A));
        Rseries = RTiOx + R0 * (1 + R0 * alphaline * pow(I_discplugserial, 2) * Rthline);
        return (Rdisc + Rplug + Rseries) * I_discplugserial;
    }

    void updateConcentration(float V_ion, float dt) {
        Nchange = (-trig / (A * lvar * 1e-9 * P_Q * zvo)) * (V_ion / 1e26);
        Nreal = Ninitreal + Nchange;
    }

    float computeVion(float V_applied, float V_schottky, float V_discplugserial) {
        if ((Nreal < Ndiscmin && V_applied > 0) | (Nreal > Ndiscmax && V_applied < 0)) { // Keep concentration Nreal in the borders of Ndiscmin and Ndiscmax
            trig = 0;
            return 0;
        } else {
            cvo = (Nplug + Nreal) / 2 * 1e26;
            if (V_applied > 0) {
                E_ion = (V_schottky + V_discplugserial * (Rdisc + Rplug) / (Rdisc + Rplug + Rseries)) / (lcell * 1e-9);
                Rtheff = Rth0 * pow(rdet/rvar, 2) * Rtheff_scaling;
                Flim = 1 - pow(Ndiscmin/Nreal, 10);
            } else {
                E_ion = V_discplugserial * Rdisc / (Rdisc + Rplug + Rseries) / (lvar * 1e-9);
                Rtheff = Rth0 * pow(rdet/rvar, 2);
                Flim = 1 - pow(Nreal/Ndiscmax, 10);
            }
            gamma = zvo * P_Q * E_ion * a / (M_PI * dWa * P_Q);
            dWamin = dWa * P_Q * (sqrt(1 - pow(gamma, 2)) - gamma * M_PI/2 + gamma * asin(gamma));
            dWamax = dWa * P_Q * (sqrt(1 - pow(gamma, 2)) + gamma * M_PI/2 + gamma * asin(gamma));
            return zvo * P_Q * cvo * a * nyo * A * (exp(-dWamin / (P_K * Treal)) - exp(dWamax / (P_K * Treal))) * Flim;
        }
    }
};
