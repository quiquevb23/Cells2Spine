#!/bin/bash

#Script to align rat reference to our fasta for Single-cell and Spatial

#We need to have samples stored in the following directories, with sample name subfolders and 
#fastqfiles subfolder in each

#single_cell_dir=/storage/gge/Quique/Cells2SpineData/Pilot/RawData/HN00224347_10X_RawData_Outs/single$
#single_cell_output=/storage/gge/Quique/Cells2SpineData/Pilot/single-cell/matrices
spatial_dir="/storage/gge/Quique/Cells2SpineData/Pilot/RawData/HN00224347_10X_RawData_Outs/spatial"
spatial_output="/storage/gge/Quique/Cells2SpineData/Pilot/spatial/matrices"
transcriptome_path=/storage/gge/Quique/rat_reference/correct_output/
echo $transcriptome_path

# Loop through the single-cell samples
for sample_dir in $spatial_dir/*; do
    sample=$(basename $sample_dir)
    echo "Processing single-cell sample: $sample"

    fastq_dir=$sample_dir/HNWKGDSXC
    image_dir=$(ls "$sample_dir/images/"*.tif 2> /dev/null)
    json_dir=$(ls "$sample_dir/json/"*.json 2> /dev/null)

    # Map sample names to areas
    case "$sample" in
        "Spatial_1")
            area="A1"
            ;;
        "Spatial_2")
            area="B1"
            ;;
        "Spatial_3")
            area="C1"
            ;;
        "Spatial_4")
            area="D1"
            ;;
        *)
            echo "Unkown sample name: $sample"
            exit 1
            ;;
    esac

    mkdir -p "$spatial_output"/"$sample"
    spaceranger count --id="$sample" \
                      --transcriptome="$transcriptome_path" \
                      --fastqs="$fastq_dir" \
                      --sample="$sample" \
                      --output-dir="$spatial_output/$sample" \
                      --create-bam=true \
                      --slide=V14J29-031 \
                      --area="$area" \
                      --image="$image_dir" \
                      --loupe-alignment="$json_dir"
done

echo "Alignment process completed."


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


