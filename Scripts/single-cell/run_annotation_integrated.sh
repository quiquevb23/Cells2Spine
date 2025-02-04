#!/bin/bash

#SBATCH --nodes=1
#SBATCH --output=annotation_int.out
#SBATCH --error=annotation_int_error.log
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=6-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --qos=medium

source /home/quiquevb/Cells2Spine/Cells2Spine/Scripts/config.sh

module load anaconda
#export PYTHONNOUSERSITE="literallyanyletters"
conda activate scanpy-env

python3 annotation_integrated.py --base_dir "$BASE_DIR" --output_base_dir "$OUTPUT_BASE_DIR"

