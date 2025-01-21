@echo off
REM Build script for Advent of Code project

REM Paths to source files
set OUTPUT=cam.exe

REM Compile command
echo Compiling...
g++ -std=c++17 -O3 ^
-I .\eigen\ ^
cam.cpp ^
main.cpp ^
-o %OUTPUT%

REM Check for success
if %errorlevel% neq 0 (
    echo Build failed.
    exit /b %errorlevel%
)

echo Build succeeded! Executable: %OUTPUT%
