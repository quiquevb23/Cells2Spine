#!/bin/bash

#Script to align rat reference to our fasta for Single-cell and Spatial

#We need to have samples stored in the following directories, with sample name subfolders and 
#fastqfiles subfolder in each

#We put now the deep seq
single_cell_dir="/storage/gge/Quique/Cells2SpineData/Pilot/RawData/HN00227613_10X_RawData_Outs"
single_cell_output="/storage/gge/Quique/Cells2SpineData/Pilot/single-cell/full_seq_join"
#spatial_dir="/storage/gge/Quique/Cells2SpineData/Pilot/spatial"
#spatial_output="$spatial_dir"/matrices
transcriptome_path=/storage/gge/Quique/rat_reference/correct_output/
echo $transcriptome_path

#The code below is for running "cellranger count" considering the 2 runs for a single matrix (correct approach)
for sample_dir in "$single_cell_dir"/*; do
    sample=$(basename "$sample_dir")
    echo "Processing single-cell sample: $sample"
    
    # Collect all fastq directories for the sample
    fastq_dirs=()
    for run_dir in "$sample_dir"/*; do
        if [ -d "$run_dir" ]; then
            fastq_dirs+=("$run_dir")
        fi
    done
    
    # Convert array to comma-separated string
    fastq_paths=$(IFS=,; echo "${fastq_dirs[*]}")
    
    # Define output directory
    output_dir="${single_cell_output}/${sample}/outs"
    mkdir -p "$output_dir"
    
    # Run Cell Ranger combining all fastq directories for the sample
    cellranger count --id="$sample" \
                     --transcriptome="$transcriptome_path" \
                     --fastqs="$fastq_paths" \
                     --sample="$sample" \
                     --create-bam=true \
                     --output-dir="$output_dir"
done


#The code below is for running "cellranger count" on each run separately
#Loop through the single-cell samples
#for sample_dir in $single_cell_dir/*; do
#    sample=$(basename $sample_dir)
#    echo "Processing single-cell sample: $sample"
#    for run_dir in $sample_dir/*; do
#        run=$(basename $run_dir)  # e.g., HNNLVDSXC or HNWKGDSXC
#        echo "Processing single-cell sample: $sample, run: $run"
        
#        output_dir="${single_cell_output}/${sample}/${run}/outs"  # Unique output directory for each run
#        mkdir -p "$output_dir"
        # Run Cell Ranger for each sequencing run
#        cellranger count --id="$sample_$run" \
#                         --transcriptome="$transcriptome_path" \
#                         --fastqs="$run_dir" \
#                         --sample="$sample" \
#                         --create-bam=true \
#                         --output-dir="$output_dir"
#    done
#done

#Loop through spatial samples
#for sample_dir in "$spatial_dir"/*; do
#    sample=$(basename "$sample_dir")
#    echo "Procecssing spatial sample: $sample"
#    mkdir -p "$spatial_output"/"$sample"
#    spaceranger count --id=$spatial_output/$sample \
#                      --transcriptome=$transcriptome_path \
#                      --fastqs=$sample_dir/fastqfiles/ \
#                      --sample=$sample \
#                      --image=$spatial_dir/images/${sample}.tiff \
#                      --slide= \ #here we need to specify slide name
#                      --area=A1 \
#done

echo "Alignment process completed."


