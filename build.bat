@echo off
REM Build script for Advent of Code project

REM Paths to source files
set OUTPUT1=RRAM_CA.exe
set OUTPUT2=RRAM_validation.exe
set OUTPUT3=RRAM_test_setup.exe

REM Compile command
echo Compiling...
g++ -std=c++17 -O3 ^
-I .\crossbar_model\eigen\ ^
-I .\crossbar_model\ ^
-I .\memristor_model\ ^
crossbar_simulator.cpp ^
nonlinear_crossbar_solver.cpp ^
main.cpp ^
.\crossbar_model\linear_crossbar_solver.cpp ^
.\memristor_model\JART_VCM_v1b_var.cpp ^
-o %OUTPUT1%

g++ -std=c++17 -O3 ^
-I .\crossbar_model\eigen\ ^
-I .\crossbar_model\ ^
-I .\memristor_model\ ^
crossbar_simulator.cpp ^
nonlinear_crossbar_solver.cpp ^
RRAM_validation.cpp ^
.\crossbar_model\linear_crossbar_solver.cpp ^
.\memristor_model\JART_VCM_v1b_var.cpp ^
-o %OUTPUT2%

g++ -std=c++17 -O3 ^
-I .\crossbar_model\eigen\ ^
-I .\crossbar_model\ ^
-I .\memristor_model\ ^
crossbar_simulator.cpp ^
nonlinear_crossbar_solver.cpp ^
RRAM_test_setup.cpp ^
.\crossbar_model\linear_crossbar_solver.cpp ^
.\memristor_model\JART_VCM_v1b_var.cpp ^
.\memristor_model\VTEAM.cpp ^
-o %OUTPUT3%

REM Check for success
if %errorlevel% neq 0 (
    echo Build failed.
    exit /b %errorlevel%
)

echo Build succeeded! Executables: %OUTPUT1%, %OUTPUT2%, and %OUTPUT3%
