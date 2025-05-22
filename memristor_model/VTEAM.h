#ifndef VTEAM_H_
#define VTEAM_H_

#include "memristor.h"

class VTEAM: public Memristor {
public:
    float vOff = 0.5;
    float vOn = -0.5;
    float kOff = 10.;
    float kOn = -100.;
    float alphaOff = 1.;
    float alphaOn = 3.;
    float Ron = 100.;
    float Roff = 1000;
    float wOn = 0;
    float wOff = 10;

    float w;
    
    void UpdateStateVariable(float v, float dt);
    float WindowFunctionOn(float w);
    float WindowFunctionOff(float w);
    float IV(float v);
    
public:
    VTEAM() {
        w = wOff;
    }

    double ApplyVoltage(double v, double dt) override;
    double GetResistance(double v) override;
    void SetWeight(bool weight) override;
};

#endif  // VTEAM_H_
