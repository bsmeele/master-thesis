@echo off
setlocal enabledelayedexpansion

:: Set the base directory where the folders are located
:: set base_dir=RRAM_validation_data\bin
set base_dir=data\bin\batch_1

:: Capture the start time in seconds (hours * 3600 + minutes * 60 + seconds)
for /f "tokens=1,2 delims=:," %%a in ("%time%") do (
    set /a start_time=%%a*3600 + %%b*60 + %%c
)

:: Print the start time (optional)
echo Start time: %time%

:: Loop through all folders with the specified pattern in the specified directory
for /d %%f in (%base_dir%\row_*_col_*) do (
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

:: Calculate the elapsed time in seconds
set /a elapsed_time=%end_time% - %start_time%

:: Handle case where the end time is earlier than the start time (overnight case)
if %elapsed_time% lss 0 (
    set /a elapsed_time+=86400  :: Add 24 hours worth of seconds (86400 seconds in a day)
)

:: Print the elapsed time in seconds
echo Elapsed time: %elapsed_time% seconds

:: Pause the script to see results before the window closes
pause
