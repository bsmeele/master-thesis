@echo off
REM Build script for Advent of Code project

REM Paths to source files
set OUTPUT=model.exe

REM Compile command
echo Compiling...
g++ -std=c++17 -g ^
JART_VCM_v1b_var.cpp ^
main.cpp ^
-o %OUTPUT%

REM Check for success
if %errorlevel% neq 0 (
    echo Build failed.
    exit /b %errorlevel%
)

echo Build succeeded! Executable: %OUTPUT%