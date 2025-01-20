@echo off
REM Build script for Advent of Code project

REM Paths to source files
set OUTPUT=RRAM_CA.exe

REM Compile command
echo Compiling...
g++ -std=c++17 -O3 ^
-I .\crossbar_model\eigen\ ^
-I .\crossbar_model\ ^
-I .\memristor_model\ ^
RRAM_crossbar_model.cpp ^
.\crossbar_model\cam.cpp ^
.\memristor_model\JART_VCM_v1b_var.cpp ^
-o %OUTPUT%

REM Check for success
if %errorlevel% neq 0 (
    echo Build failed.
    exit /b %errorlevel%
)

echo Build succeeded! Executable: %OUTPUT%
