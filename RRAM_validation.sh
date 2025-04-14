#!/bin/bash

# Set the base directory where the folders are located
# base_dir="./weights_inputs/bin"
base_dir="./data/bin/batch_1"

# Capture the start ime in seconds (hours * 3600 + minutes * 60 + seconds)
start_time=$(date +%s)

# Print start time
echo "Start time: $(date)"

# Loop through all folders with the specified pattern in the specified directory
for folder in "$base_dir"/row_*_col_*_*; do
    echo "Processing folder: $folder"

    # Call your C++ program with the folder as an argument
    ./RRAM_validation.out "$folder"
done

# Capture the end time in seconds
end_time=$(date +%s)

# Print the end time
echo "End time: $(date)"
