#include "VTEAM.h"

#include <cmath>
#include <iostream>
#include <sstream>

void VTEAM::UpdateStateVariable(float v, float dt) {
    float dw = 0;
    if (vOff <= v) {
        dw = kOff * pow(v/vOff - 1, alphaOff) * WindowFunctionOff(w);
    } else if (v <= vOn) {
        dw = kOn * pow(v/vOn - 1, alphaOn) * WindowFunctionOn(w);
    }

    w += dw * dt;
}

float VTEAM::WindowFunctionOff(float w) {
    if (w < wOff) { return 1.; }
    else { return 0.; }
}

float VTEAM::WindowFunctionOn(float w) {
    if (w > wOn) { return 1.; }
    else { return 0.; }
}

float VTEAM::IV(float v) {
    return pow(Ron + (Roff - Ron)/(wOff - wOn) * (w - wOn), -1) * v;
}

double VTEAM::ApplyVoltage(double v, double dt) {
    UpdateStateVariable(v, dt);
    return IV(v);
}

double VTEAM::GetResistance(double v) {
    return Ron + (Roff - Ron)/(wOff - wOn) * (w - wOn);
}

void VTEAM::SetWeight(bool weight) {
    if (weight) { w = wOn; }
    else { w = wOff; }
}

std::string VTEAM::GetParams() {
    std::ostringstream params;
    params << "Model: VTEAM" << std::endl
        << "vOff:" << vOff << std::endl
        << "vOn: " << vOn << std::endl
        << "kOff: " << kOff << std::endl
        << "kOn: " << kOn << std::endl
        << "alphaOff: " << alphaOff << std::endl
        << "alphaOn: " << alphaOn << std::endl
        << "Ron: " << Ron << std::endl
        << "Roff: " << Roff << std::endl
        << "wOn: " << wOn << std::endl
        << "wOff: " << wOff << std::endl
        << "d: " << d << std::endl;
    return params.str();
}

