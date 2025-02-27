#!/bin/bash
# Build script for master thesis code


# Paths to source files
OUTPUT="cam.out"

# Compile command
echo "Compiling..."
g++ -std=c++17 -g \
-I eigen/ \
cam.cpp \
main.cpp \
-o $OUTPUT

# Check for success
if [ $? -ne 0 ]; then
    echo "Build failed."
    exit 1
fi

echo "Build succeeded! Executable: $OUTPUT"
