#ifndef DUMMY_MEMRISTOR_H_
#define DUMMY_MEMRISTOR_H_

#include "memristor.h"

class Dummy: public Memristor {
public:
    const float Ron = 100;
    const float Roff = 10000;
    float R;

    Dummy() {  R = Roff; }

    double ApplyVoltage(double V_applied, double dt) override {
        return V_applied / R;
    };
    double GetResistance(double V_applied) {
        return R;
    };
    void SetWeight(bool weight) {
        if (weight) { R = Ron; }
        else { R = Roff; }
    };
};

#endif  // DUMMY_MEMRISTOR_H_
