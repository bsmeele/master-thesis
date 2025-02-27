#!/bin/bash
# Build script for master thesis code


# Paths to source files
OUTPUT1="RRAM_CA.out"
OUTPUT2="RRAM_validation.out"

# Compile command
echo "Compiling..."
g++ -std=c++17 -O3 \
-I crossbar_model/eigen/ \
-I crossbar_model/ \
-I memristor_model/ \
crossbar_simulator.cpp \
nonlinear_crossbar_solver.cpp \
main.cpp \
crossbar_model/linear_crossbar_solver.cpp \
memristor_model/JART_VCM_v1b_var.cpp \
-o $OUTPUT1

g++ -std=c++17 -O3 \
-I crossbar_model/eigen/ \
-I crossbar_model/ \
-I memristor_model/ \
crossbar_simulator.cpp \
nonlinear_crossbar_solver.cpp \
RRAM_validation.cpp \
crossbar_model/linear_crossbar_solver.cpp \
memristor_model/JART_VCM_v1b_var.cpp \
-o $OUTPUT2

# Check for success
if [ $? -ne 0 ]; then
    echo "Build failed."
    exit 1
fi

echo "Build succeeded! Executable: $OUTPUT1 and $OUTPUT2"
