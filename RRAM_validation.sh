#!/bin/bash

# Set the base directory where the folders are located
base_dir="./weights_inputs/bi"

# Capture the start ime in seconds (hours * 3600 + minutes * 60 + seconds)
start_time=$(date +%s)

# Print start time
echo "Start time: $(date)"

# Loop through all folders with the specified pattern in the specified directory
for folder in "$base_dir"/row_*_col_*_*; do
    echo "Processing folder: $folder"

    # Call your C++ program with the folder as an argument
    ./RRAM_validation "$folder"
done

# Capture the end time in seconds
end_time=$(date +%s)

# Print the end time
echo "End time: $(date)"

# Calculate the elapsed time in seconds
elapsed_time=$((end_time - start_time))

# Handle case where the end time is earlier than the start time (overnight case)
if [ $elapsed_time -lt 0 ]; then
    elapsed_time=$((elapsed_time + 86400))  # Add 24 hours worth of seconds (86400 seconds in a day)
fi

# Print the elapsed time in seconds
echo "Elapsed time: $elapsed_time seconds"
