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

// TODO:
//   Variability model
//   Add changing fitting parameters
//   Fix early reset
//     I'm still working out why this happens. It seems to be caused by a sharp increase in I_ion, which is caused by Treal, which is caused by V_schottky/Ischottky, which might be caused by phibn
//     The point where it changes seems to also be when (V_schottky < phibn0 - phin) is triggered to be false
//     I'm yet unsure if this causes the bug or is a result of the bug
//     However, it should be noted that phibn for V_schottky=0.0799 is 0.0621074, while phibn for V_schottky = 0.08 (which triggers it condition) is 0.18, which is phibn0
//     This makes me believe that this discontinuity is the culprit of the early state transition, however I'm unsure if the discontinuity is intentional
//   Fix jitter during reset
//     This bug seems to have the same origin as the early reset

// Constants
#ifndef M_PI
#define M_PI 3.1415927  // Define M_PI if not defined
#endif
const double P_Q = 1.6022e-19;    // Elementary charge [C]
const double P_K = 1.38065e-23;   // Boltzman constant [J/K]
const double P_EPS0 = 8.6549e-12; // Permittivity of a vacuum [F/m]
const double P_H = 6.626e-34;     // Planck constant [Js]

class JART_VCM_v1b_var {
    private:
    public:
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

        double Nreal; // oxygen vacancy concentration of the disc region [nm]

        void updateFilamentArea() {
            A = M_PI * pow(rvar, 2);
        }

        void updateTemperature(double V_schottky, double V_discplugserial, double I_schottky) {
            // if (V_prev > 0.63 && V_prev < 0.642) {
            //     std::cout << Treal << ' ' << I_schottky << ' ' << Rtheff << ' ' << V_schottky << ' ' << V_discplugserial << ' ' << Rdisc << ' ' << Rplug << ' ' << Rseries << ' ' << Nreal << std::endl;
            // }
            Treal = I_schottky * (V_schottky + V_discplugserial * (Rdisc + Rplug) / (Rdisc + Rplug + Rseries)) * Rtheff + T0;
        }

        double computeSchottkyCurrent(double V_schottky, bool print=false) {
            double phibn;

            if (V_schottky < phibn0 - phin) {
                double psi = phibn0 - phin - V_schottky;
                phibn = phibn0 - sqrt(sqrt(pow(P_Q, 3) * zvo * Nreal * 1e26 * psi / (8 * pow(M_PI, 2) * (pow(epsphib_eff, 3)))));
                if (phibn < 0) { phibn = 0; }
            } else { phibn = phibn0; }

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

        void updateResistance(double I_discplugserial) {
            Rdisc = lvar * 1e-9 / (Nreal * 1e26 * zvo * P_Q * un * A);
            Rplug = ((lcell - lvar) * 1e-9 / (Nplug * 1e26 * zvo * P_Q * un * A));
            Rseries = RTiOx + R0 * (1 + R0 * alphaline * pow(I_discplugserial, 2) * Rthline);
        }

        void updateConcentration(double I_ion, double dt) {
            double Nchange = (-trig / (A * lvar * 1e-9 * P_Q * zvo) * I_ion / 1e26) * dt;
            // Nreal = Ninitreal + Nchange;
            Nreal += Nchange;
            if (Nreal > Ndiscmax) { Nreal = Ndiscmax; }
            else if (Nreal < Ndiscmin) { Nreal = Ndiscmin; }
        }

        double computeIonCurrent(double V_applied, double V_schottky, double V_discplugserial) {
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

    // public:
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
            updateResistance(0);
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

            double V_low = 0;
            double V_high = V_applied;

            double V_schottky;
            double I_schottky;
            double V_discplugserial;

            double tresh = 0.64;
            double tresh2 = 0.72;
            int i = 0;
            while (true) {
                V_schottky = (V_low + V_high) / 2;
                if (i > 1e3) {
                    std::cout << "Iteration limit reached" << std::endl;
                    int test = 1/0;
                }
                // if (V_applied > tresh) {
                //     std::cout << "V applied: " << V_applied << std::endl;
                //     std::cout << "V schottky: " << V_schottky << std::endl;
                // }
                if (std::isinf(V_schottky) || std::isinf(I_schottky) || std::isinf(V_discplugserial)) {
                    std::cout << "inf detected" << std::endl;
                    int test = 1 / 0;
                }
                if (std::isnan(V_schottky) || std::isnan(I_schottky) || std::isnan(V_discplugserial)) {
                    std::cout << "nan detected" << std::endl;
                    int test = 1 / 0;
                }

                I_schottky = computeSchottkyCurrent(V_schottky);

                // if (V_applied > tresh) {
                //     std::cout << "I schottky: " << I_schottky << std::endl;
                // }

                updateResistance(I_schottky);

                V_discplugserial = (Rdisc + Rplug + Rseries) * I_schottky;

                // if (V_applied > tresh) {
                //     std::cout << "V discplugserial: " << V_discplugserial << std::endl;
                //     std::cout << "Criterion: " << fabs((V_applied - V_discplugserial) - V_schottky) << std::endl;
                // }

                double err = V_applied - V_discplugserial - V_schottky;
                if (fabs(err) < 1e-6) { break; }
                if (V_applied < 0) {
                    if (err < 0) { V_low = V_schottky; }
                    else { V_high = V_schottky; }
                } else {
                    if (err > 0) { V_low = V_schottky; }
                    else { V_high = V_schottky; }
                }

                // if (V_applied > tresh) {
                //     std::cout << "V_applied - V_discplugserial: " << V_applied - V_discplugserial << std::endl;
                // }

                i += 1;
                
                // if (V_applied > tresh) {
                //     std::cout << std::endl;
                // }
            }

            if (V_prev > 0.157 && V_prev < 0.162) {
                double tmp = computeSchottkyCurrent(V_schottky, true);
                // std::cout << V_applied << ' ' << V_schottky << ' ' << tmp << ' ' << Nreal << std::endl;
            }

            double I_ion = computeIonCurrent(V_applied, V_schottky, V_discplugserial);
            // if (V_applied > tresh && V_applied < tresh2) {
            //     std::cout << "dt: " << dt << ", V_applied: " << V_applied << ", I_ion: " << I_ion << ", Nreal: " << Nreal << std::endl;
            // }
            double N_before = Nreal;
            updateConcentration(I_ion, dt);

            // Force smaller time steps during abrupt switching
            if (fabs(Nreal - N_before) > 1e-1) {
                Nreal = N_before;
                for (int i = 0; i < 10; i++) {
                    apply_voltage(V_applied, dt/10.);
                }
            } else {
                updateTemperature(V_schottky, V_discplugserial, I_schottky);
            }

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

    outFile << "t V I Nreal Treal Rschottky Rdisc Rplug Rseries" << std::endl;

    for (int i = 1; i < V_wave.size(); i++) {
        // std::cout << i << std::endl;
        double dv = (V_wave[i][0] - V) / ((V_wave[i][1] - t) / dt);
        while (t < V_wave[i][1]) {
            double I;
            I = memristor.apply_voltage(V, dt);
            // std::cout << t << "s" << ": " << V << " V" << ", " << I << " A, " << "N: " << memristor.Nreal << std::endl;
            if (std::isnan(I)) { return 1; }
            outFile << t << " " << V << " " << I << " " << memristor.Nreal << " " << memristor.Treal << " " << (V - (memristor.Rdisc + memristor.Rplug + memristor.Rseries) * I)/I << " " << memristor.Rdisc << " " << memristor.Rplug << " " << memristor.Rseries << std::endl;
            V += dv;
            t += dt;
        }
    }

    outFile.close();

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
    //     memristor.computeSchottkyCurrent(v, true);
    //     // std::cout << v << ' ' << memristor.computeSchottkyCurrent(v) << std::endl;
    // }
}
