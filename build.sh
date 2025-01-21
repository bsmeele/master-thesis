#!/bin/bash
# Build script for master thesis code


# Paths to source files
OUTPUT="RRAM_CA"

# Compile command
echo "Compiling..."
g++ -std=c++17 -g \
-I crossbar_model/eigen/ \
-I crossbar_model/ \
-I memristor_model/ \
RRAM_crossbar_model.cpp \
crossbar_model/cam.cpp \
memristor_model/JART_VCM_v1b_var.cpp \
-o $OUTPUT

# Check for success
if [ $? -ne 0 ]; then
    echo "Build failed."
    exit 1
fi

echo "Build succeeded! Executable: $OUTPUT"
