#!/bin/bash

#Script to run CellBender: it requires > 24h, need to be run on individual samples and unfiltered data

output_dir="/storage/gge/Quique/Cells2SpineData/Pilot/single-cell/full_seq"

# List of input files to process
input_files=(
    "$output_dir/Single_Cell_1/HWYVKDSXC/outs/outs/raw_feature_bc_matrix.h5"
)

#Loop through the single-cell samples
for input_file in "${input_files[@]}"; do
    # Extract the sample name from the path
    sample=$(basename "$(dirname "$(dirname "$input_file")")") # Gets the sample folder name (e.g., Single_Cell_1)
    echo "Processing single-cell sample: $sample"
    
    # Define the output file path
    output_file="${input_file/raw_feature_bc_matrix.h5/cellbender_filtered.h5}"
    
    # Run CellBender
    cellbender remove-background \
        --input "$input_file" \
        --output "$output_file"
done


