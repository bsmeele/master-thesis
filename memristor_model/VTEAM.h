#ifndef VTEAM_H_
#define VTEAM_H_

#include "memristor.h"

class VTEAM: public Memristor {
public:
    float vOff = 0.02;
    float vOn = -0.2;
    float kOff = 5e-4;
    float kOn = -10.;
    float alphaOff = 1.;
    float alphaOn = 3.;
    float Ron = 50.;
    float Roff = 1000;
    float wOn = 0;
    float wOff = 3e-9;
    float d = 3e-9;

    float w;
    
    void UpdateStateVariable(float v, float dt);
    float WindowFunctionOn(float w);
    float WindowFunctionOff(float w);
    float IV(float v);
    
public:
    VTEAM() {
        w = wOn;
    }

    double ApplyVoltage(double v, double dt) override;
    double GetResistance(double v) override;
    void SetWeight(bool weight) override;
};

#endif  // VTEAM_H_
