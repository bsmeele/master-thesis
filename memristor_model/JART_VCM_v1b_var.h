#ifndef JART_VCM_v1b_var_H_
#define JART_VCM_v1b_var_H_

#include "memristor.h"

#include <array>
#include <vector>
#include <iostream>

#ifndef M_PI
#define M_PI 3.1415927  // Define M_PI if not defined
#endif
const double P_Q = 1.6022e-19;    // Elementary charge [C]
const double P_K = 1.38065e-23;   // Boltzman constant [J/K]
const double P_EPS0 = 8.6549e-12; // Permittivity of a vacuum [F/m]
const double P_H = 6.626e-34;     // Planck constant [Js]

class JART_VCM_v1b_var: public Memristor {
    private:
    public:
        // ----- Pyisical constants do not change! -----
        const double Arichardson = 6.01e5; // Richardson's constant [A/m^2K^2]
        const double mdiel = 9.10938e-31;  // electron rest mass [kg]
        const double zvo = 2;              // oxygen vacancy charge number
        const double T0 = 293;             // ambient temperature [K]

        // ----- Fitting parameters -----
        const double eps = 17;              // from [10:25], static hafnium oxide permittivity
        const double epsphib = 5.5;         // hafnium oxide permittivity related to image force barrier lowering
        // const double phibn0 = 0.18;         // from [0.1:0.5], nominal schottky barrier height [eV]
        const double phibn0 = 0.18;
        const double phin = 0.1;            // from [0.1:0.3], energy level difference between the Fermi level in the oxide and the oxide conduction band edge [eV]
        const double un = 4e-6;             // from [1e-6:1e-5], electron mobility [m^2/Vs]
        const double Ndiscmax = 20.;         // from [0.001:1100], maximum oxygen vacancy concentration in the disc [10^26/m^3]
        // const double Ndiscmin = 0.008;      // from [0.0001:100], minimum oxygen vacancy concentration in the disc [10^26/m^3]
        // const double Ninit = 0.008;         // from [0.0001:1000], initial oxygen vacancy concentration in the disc [10^26/m^3]
        const double Ndiscmin = 0.0001;      // from [0.0001:100], minimum oxygen vacancy concentration in the disc [10^26/m^3]
        const double Ninit = 0.0001;         // from [0.0001:1000], initial oxygen vacancy concentration in the disc [10^26/m^3]
        // const double Ninit = 20;
        // const double Ndiscmin = 0.0001;
        // const double Ninit = 20.;
        const double Nplug = 20;            // from [0.001:100], oxygen vacancy concentration in the plug [10^26/m^3]
        const double a = 0.25e-9;           // from [0.1e-9:1e-9], ion hopping distance [m]
        const double nyo = 2e13;            // from [1e10:1e14], attempt frequency [Hz]
        const double dWa = 1.35;            // from [0.8:1.5], activation energy [eV]
        // const double Rth0 = 15.72e6;        // from [1e6:20e6], thermal resistance of the Hafnium Oxide [K/W]
        const double Rth0 = 1e7;
        const double rdet = 45e-9;          // from [5e-9:100e-9], radius of the filament [m]
        const double rnew = 45e-9;          // from [5e-9:100e-9], radius of the filament [m]
        const double lcell = 3;             // from [2:5], length of disc and plug region [nm]
        const double ldet = 0.4;            // from [0.1:5], length of the disc region [nm]
        const double lnew = 0.4;            // from [0.1:5], length of the disc region [nm]
        const double Rtheff_scaling = 0.27; // from [0.1:1], scaling factor for RESET
        const double RTiOx = 650;           // from [0:5000], series resistance of the TiOx layer [Ohm]
        const double R0 = 719.2437;         // Resistance at T0 [Ohm]
        const double Rthline = 90471.47;    // thermal conductivity of the Platinum and Titanium [W/mK]
        const double alphaline = 3.92e-3;   // temperature coefficient [1/K]

        // const double RTiOx = 10;
        // const double Rthline = 1000000;
        // const double alphaline = 1e-1;

        double eps_eff;     // static hafnium oxide permittivity
        double epsphib_eff; // hafnium oxide permittivity related to image force barrier lowering

        // ----- State variable -----

        // ----- Other internal variables -----
        int trig;         // Used to signify certain voltage crossings and limit the state variable
        double Ninitreal; // Not sure what this does, it is set to Ninit on startup and never changed again
        double rvar;      // Radius of the fillament used in calculations, updated with the variability model
        double rold;      // Radius of the fillament, only used in the variability model to update rvar
        double lvar;      // Length of the disc region used in calculations, updated with the variability model
        double lold;      // Length of the disc region, only used in the variability model to update lvar
        double Nold;      // Oxygen vacancy of the disc region, only used in the variability model
        double Treal;     // Homogeneous filament temperature [K]
        double A;         // Cross section of the filament area
        double Rdisc;     // Resistance of the disc region
        double Rplug;     // Resistance of the plug region
        double Rseries;   // Resistance of the series section
        double Rtheff;    // Thermal resistance
        double V_prev;    // |revious applied voltage, used to check for voltage crossings
        double V_schottky_prev;  // Schottky voltage calculated during the previous iteration

        double Nreal; // oxygen vacancy concentration of the disc region [nm]

        double tmp;

        void UpdateFilamentArea();
        void UpdateTemperature(double V_schottky, double V_discplugserial, double I_schottky);
        double ComputeSchottkyCurrent(double V_schottky);
        void UpdateResistance(double I_discplugserial);
        void UpdateConcentration(double I_ion, double dt);
        double ComputeIonCurrent(double V_applied, double V_schottky, double V_discplugserial);
        std::array<double, 3> SolveFixedpoint(double Vguess, double V_applied);

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
            V_prev = 0;
            V_schottky_prev = 0;
            UpdateFilamentArea();
            UpdateResistance(0);
            UpdateTemperature(0, 0, 0);

            tmp = 0;
        }
        double ApplyVoltage(double V_applied, double dt) override;
        double GetResistance(double V_applied) override;
        void SetWeight(bool weight) override;
        std::string GetParams() override;
};

#endif  // JART_VCM_v1b_var_H_
