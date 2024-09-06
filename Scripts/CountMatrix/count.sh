#!/bin/bash

#Script to align rat reference to our fasta for Single-cell and Spatial

#We need to have samples stored in the following directories, with sample name subfolders and 
#fastqfiles subfolder in each

single_cell_dir="/storage/gge/Quique/Cells2SpineData/Pilot/single-cell"
single_cell_output="$single_cell_dir"/matrices
spatial_dir="/storage/gge/Quique/Cells2SpineData/Pilot/spatial"
spatial_output="$spatial_dir"/matrices
transcriptome_path="/storage/gge/Quique/rat_reference/correct_output/"

mkdir -p "$single_cell_output"
mkdir -p "$spatial_output"

#Loop through the single-cell samples
for sample_dir in "$single_cell_dir"/*; do
    sample=$(basename "$sample_dir")
    echo "Processing single-cell sample: $sample"
    cellranger count --id=$single_cell_output/$sample \
                     --transcriptome=$transcriptome_path \
                     --fastqs=$sample_dir/fastqfiles/ \
                     --sample=$sample \
                     --create-bam true
done


#Loop through spatial samples
for sample_dir in "$spatial_dir"/*; do
    sample=$(basename "$sample_dir")
    echo "Procecssing spatial sample: $sample"
    mkdir -p "$spatial_output"/"$sample"
    spaceranger count --id=$spatial_output/$sample \
                      --transcriptome=$transcriptome_path \
                      --fastqs=$sample_dir/fastqfiles/ \
                      --sample=$sample \
                      --image=$spatial_dir/images/${sample}.tiff \
                      --slide= \ #here we need to specify slide name
                      --area=A1 \
done

echo "Alignment process completed."
