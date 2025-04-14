#!/bin/bash
# Build script for master thesis code


# Paths to source files
OUTPUT="model.out"

# Compile command
echo "Compiling..."
g++ -std=c++17 -g \
JART_VCM_v1b_var.cpp \
main.cpp \
-o $OUTPUT

# Check for success
if [ $? -ne 0 ]; then
    echo "Build failed."
    exit 1
fi

echo "Build succeeded! Executable: $OUTPUT"
