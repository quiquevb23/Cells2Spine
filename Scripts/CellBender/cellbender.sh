#!/bin/bash

#Script to run CellBender: it requires > 24h, need to be run on individual samples and unfiltered data

single_cell_dir="/storage/gge/Quique/Cells2SpineData/Pilot/single-cell"
single_cell_output="$single_cell_dir"/matrices
spatial_dir="/storage/gge/Quique/Cells2SpineData/Pilot/spatial"
spatial_output="$spatial_dir"/matrices
transcriptome_path="/storage/gge/Quique/rat_reference/correct_output/"

mkdir -p "$single_cell_output"
mkdir -p "$spatial_output"

#Loop through the single-cell samples
for sample_dir in "$single_cell_output"/*; do
    sample=$(basename "$sample_dir") #retain only first element before underscore
    echo "Processing single-cell sample: $sample"
    cellbender remove-background \
            --input $sample_dir/raw*
            --output $sample_dir/
done

