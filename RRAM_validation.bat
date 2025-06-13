@echo off
setlocal enabledelayedexpansion

:: Set the base directory where the folders are located
:: set base_dir=RRAM_validation_data\bin
@REM set base_dir=data\bin\batch_1
set base_dir=RRAM_validation_data\64x64\save_c_arrays

:: Capture the start time in seconds (hours * 3600 + minutes * 60 + seconds)
for /f "tokens=1,2 delims=:," %%a in ("%time%") do (
    set /a start_time=%%a*3600 + %%b*60 + %%c
)

:: Print the start time (optional)
echo Start time: %time%

:: Loop through all folders with the specified pattern in the specified directory
@REM for /d %%f in (%base_dir%\row_*_col_*) do (
for /d %%f in (%base_dir%\batch_*) do (
    set folder=%%f
    echo Processing folder: !folder!
    
    :: Call your C++ program with the folder as an argument
    .\RRAM_validation.exe "!folder!"
)

:: Capture the end time in seconds (hours * 3600 + minutes * 60 + seconds)
for /f "tokens=1,2 delims=:," %%a in ("%time%") do (
    set /a end_time=%%a*3600 + %%b*60 + %%c
)

:: Print the start time (optional)
echo Start time: %time%

:: Print the end time (optional)
echo End time: %time%

:: Pause the script to see results before the window closes
pause
