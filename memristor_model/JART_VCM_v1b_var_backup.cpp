#include <cmath>
#include <vector>
#include <array>
#include <fstream>
#include <iostream>

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



// Constants
#ifndef M_PI
#define M_PI 3.1415927  // Define M_PI if not defined
#endif
const double P_Q = 1.6022e-19;    // Elementary charge [C]
const double P_K = 1.38065e-23;   // Boltzman constant [J/K]
const double P_EPS0 = 8.6549e-12; // Permittivity of a vacuum [F/m]
const double P_H = 6.626e-34;     // Planck constant [Js]

// TODO:
//   Variability model
//   Force small time steps during abrupt switching?
//   Add changing fitting parameters

class JART_VCM_v1b_var {
    private:
        // ----- Pyisical constants do not change! -----
        const double Arichardson = 6.01e5; // Richardson's constant [A/m^2K^2]
        const double mdiel = 9.10938e-31;  // electron rest mass [kg]
        const double zvo = 2;              // oxygen vacancy charge number
        const double T0 = 293;             // ambient temperature [K]

        // ----- Fitting parameters -----
        double eps = 17;              // from [10:25], static hafnium oxide permittivity
        double epsphib = 5.5;         // hafnium oxide permittivity related to image force barrier lowering
        double phibn0 = 0.18;         // from [0.1:0.5], nominal schottky barrier height [eV]
        double phin = 0.1;            // from [0.1:0.3], energy level difference between the Fermi level in the oxide and the oxide conduction band edge [eV]
        double un = 4e-6;             // from [1e-6:1e-5], electron mobility [m^2/Vs]
        double Ndiscmax = 20;         // from [0.001:1100], maximum oxygen vacancy concentration in the disc [10^26/m^3]
        double Ndiscmin = 0.008;      // from [0.0001:100], minimum oxygen vacancy concentration in the disc [10^26/m^3]
        double Ninit = 0.008;         // from [0.0001:1000], initial oxygen vacancy concentration in the disc [10^26/m^3]
        double Nplug = 20;            // from [0.001:100], oxygen vacancy concentration in the plug [10^26/m^3]
        double a = 0.25e-9;           // from [0.1e-9:1e-9], ion hopping distance [m]
        double nyo = 2e13;            // from [1e10:1e14], attempt frequency [Hz]
        double dWa = 1.35;            // from [0.8:1.5], activation energy [eV]
        double Rth0 = 15.72e6;        // from [1e6:20e6], thermal resistance of the Hafnium Oxide [K/W]
        double rdet = 45e-9;          // from [5e-9:100e-9], radius of the filament [m]
        double rnew = 45e-9;          // from [5e-9:100e-9], radius of the filament [m]
        double lcell = 3;             // from [2:5], length of disc and plug region [nm]
        double ldet = 0.4;            // from [0.1:5], length of the disc region [nm]
        double lnew = 0.4;            // from [0.1:5], length of the disc region [nm]
        double Rtheff_scaling = 0.27; // from [0.1:1], scaling factor for RESET
        double RTiOx = 650;           // from [0:5000], series resistance of the TiOx layer [Ohm]
        double R0 = 719.2437;         // Resistance at T0 [Ohm]
        double Rthline = 90471.47;    // thermal conductivity of the Platinum and Titanium [W/mK]
        double alphaline = 3.92e-3;   // temperature coefficient [1/K]

        double eps_eff;     // static hafnium oxide permittivity
        double epsphib_eff; // hafnium oxide permittivity related to image force barrier lowering

        // ----- State variable -----
        double Nreal; // oxygen vacancy concentration of the disc region [nm]

        // ----- Other internal variables -----
        int trig;        // Used to signify certain voltage crossings and limit the state variable
        double Ninitreal; // Not sure what this does, it is set to Ninit on startup and never changed again
        double rvar;      // radius of the fillament used in calculations, updated with the variability model
        double rold;      // radius of the fillament, only used in the variability model to update rvar
        double lvar;      // length of the disc region used in calculations, updated with the variability model
        double lold;      // length of the disc region, only used in the variability model to update lvar
        double Nold;      // oxygen vacancy of the disc region, only used in the variability model
        double Treal;     // homogeneous filament temperature [K]
        double A;         // cross section of the filament area
        double Rdisc;     // resistance of the disc region
        double Rplug;     // resistance of the plug region
        double Rseries;   // resistance of the series section
        double Rtheff;    // thermal resistance
        double V_prev;    // previous applied voltage, used to check for voltage crossings

        void setFittingParameters() {

        }

        void updateFilamentArea() {
            A = M_PI * pow(rvar, 2);
        }

        void updateTemperature(double V_schottky, double V_discplugserial, double I_schottky) {
            Treal = I_schottky * (V_schottky + V_discplugserial * (Rdisc + Rplug) / (Rdisc + Rplug + Rseries)) * Rtheff + T0;
        }

        double computeSchottkyCurrent(double V_schottky) {
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
                return -A * ((Arichardson * Treal) / P_K) * sqrt(M_PI * W00 * P_Q * (fabs(V_schottky) + phibn / pow(cosh(W00/(P_K * Treal)),2)))
                * exp(-P_Q * phibn / W0) * (exp(P_Q * fabs(V_schottky) / epsprime) - 1);
            } else { // Schottky TE RESET direction
                return A * Arichardson * pow(Treal, 2) * exp(-phibn * P_Q / (P_K * Treal)) * (exp(P_Q / (P_K * Treal) * V_schottky) - 1);
            }
        }

        void updateResistance(double I_discplugserial) {
            Rdisc = lvar * 1e-9 / (Nreal * 1e26 * zvo * P_Q * un * A);
            Rplug = ((lcell - lvar) * 1e-9 / (Nplug * 1e26 * zvo * P_Q * un * A));
            Rseries = RTiOx = R0 * (1 + R0 * alphaline * pow(I_discplugserial, 2) * Rthline);
        }

        void updateConcentration(double I_ion, double dt) {
            double Nchange = (-trig / (A * lvar * 1e-9 * P_Q * zvo) * I_ion / 1e26) * dt;
            Nreal = Ninitreal + Nchange;
        }

        double computeIonCurrent(double V_applied, double V_schottky, double V_discplugserial) {
            if ((Nreal < Ndiscmin && V_applied > 0) | (Nreal > Ndiscmax && V_applied < 0)) { // Keep concentration Nreal in the borders of Ndiscmin and Ndiscmax
                trig = 0;
                return 0;
            } else {
                double cvo = (Nplug + Nreal) / 2 * 1e26;
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
                return zvo * P_Q * cvo * a * nyo * A * (exp(-dWamin / (P_K * Treal)) - exp(dWamax / (P_K * Treal))) * Flim;
            }
        }

    public:
        JART_VCM_v1b_var() {
            Nreal = Ninit;
            Ninitreal = Ninit;
            trig = 1;
            rvar = rnew;
            rold = rnew;
            lvar = lnew;
            lold = lnew;
            Nold = Ninit;
            eps_eff = eps * P_EPS0;
            epsphib_eff = epsphib * P_EPS0;
            updateFilamentArea();
            updateTemperature(0, 0, 0);
        }

        double apply_voltage(double V_applied, double dt) {
            // Check voltage crossings
            if ((V_prev > -1.5e-5 && V_applied <= -1.5e-5) || (V_prev < 1.5e-5 && V_applied >= 1.5e-5)) {
                rold = rvar;
                lold = lvar;
                Nold = Nreal;
                trig = 1;
            }
            V_prev = V_applied;

            double V_schottky = V_applied;
            double I_schottky;
            double V_discplugserial;

            while (true) {
                std::cout << "V schottky: " << V_schottky << std::endl;
                if (std::isinf(V_schottky)) { int test = 1 / 0; }

                I_schottky = computeSchottkyCurrent(V_schottky);

                std::cout << "I schottky: " << I_schottky << std::endl;

                updateResistance(I_schottky);

                V_discplugserial = (Rdisc + Rplug + Rseries) * I_schottky;

                std::cout << "V discplugserial: " << V_discplugserial << std::endl;

                double criterion = fabs((V_applied - V_discplugserial) - V_schottky);
                if (criterion < 1e-6) {
                    V_schottky = V_applied - V_discplugserial;
                    break;
                }
                V_schottky = V_applied - V_discplugserial;
            }

            double I_ion = computeIonCurrent(V_applied, V_schottky, V_discplugserial);
            updateConcentration(I_ion, dt);
            updateTemperature(V_schottky, V_discplugserial, I_schottky);

            return I_schottky;
        }
};

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

    // std::ofstream outFile("Iout.txt");

    // if (!outFile) {
    //     std::cout << "No out file" << std::endl;
    //     return 1;
    // }

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

    for (int i = 1; i < V_wave.size(); i++) {
        double dv = (V_wave[i][0] - V) / ((V_wave[i][1] - t) / dt);
        while (t < V_wave[i][1]) {
            double I;
            // std::cout << t << ": " << V << " V";
            I = memristor.apply_voltage(V, dt);
            // std::cout << ", " << I << " A" << std::endl;
            if (std::isnan(I)) { return 1; }
            // outFile << I << std::endl;
            V += dv;
            t += dt;
        }
    }


    // outFile.close();
}