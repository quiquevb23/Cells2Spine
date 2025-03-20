#!/bin/bash

#SBATCH --nodes=1
#SBATCH --output=compare_ann.out
#SBATCH --error=compare_ann_error.log
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --time=6-00:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --qos=medium

#source /home/quiquevb/Cells2Spine/Cells2Spine/Scripts/config.sh
source /home/quiquevb/Cells2Spine/Cells2Spine/Scripts/config_cellbender.sh

module load anaconda
#export PYTHONNOUSERSITE="literallyanyletters"
conda activate scanpy-env

python3 compare_ann.py --base_dir "$BASE_DIR" --output_base_dir "$OUTPUT_BASE_DIR"
