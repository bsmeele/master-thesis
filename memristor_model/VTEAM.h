#ifndef VTEAM_H_
#define VTEAM_H_

#include "memristor.h"
#include "string"

class VTEAM: public Memristor {
public:
    const float vOff = 0.02;
    const float vOn = -0.2;
    const float kOff = 5e-4;
    const float kOn = -10.;
    const float alphaOff = 1.;
    const float alphaOn = 3.;
    const float Ron = 50.;
    const float Roff = 1000;
    const float wOn = 0;
    const float wOff = 3e-9;
    const float d = 3e-9;

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
    std::string GetParams() override;
};

#endif  // VTEAM_H_
