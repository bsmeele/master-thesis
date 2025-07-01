#ifndef DUMMY_NONLINEAR_MEMRISTOR_H_
#define DUMMY_NONLINEAR_MEMRISTOR_H_

#include "memristor.h"
#include <sstream>
#include <cmath>

class DummyNonlinear: public Memristor {
public:
    const float Ron = 1000;
    const float Roff = 100000;
    const float a = 1.;
    float R;

    DummyNonlinear() {  R = Roff; }

    double ApplyVoltage(double V_applied, double dt) override {
        return V_applied / ((1 + a*std::fabs(V_applied))*R);
    };
    double GetResistance(double V_applied) {
        return (1 + a*std::fabs(V_applied))*R;
    };
    void SetWeight(bool weight) {
        if (weight) { R = Ron; }
        else { R = Roff; }
    };
    std::string GetParams() override {
        std::ostringstream params;
        params << "Model: Dummy Nonlinear" << std::endl
            << "Ron: " << Ron << std::endl
            << "Roff: " << Roff << std::endl
            << "a: " << a << std::endl;
        return params.str();
    }
};

#endif  // DUMMY_NONLINEAR_MEMRISTOR_H_
