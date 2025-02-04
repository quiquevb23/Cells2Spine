#!/bin/bash

#SBATCH --nodes=1
#SBATCH --output=xgboost_annot.out
#SBATCH --error=xgboost_annot_error.log
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --time=6-00:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --qos=medium

source /home/quiquevb/Cells2Spine/Cells2Spine/Scripts/config.sh

module load anaconda
#export PYTHONNOUSERSITE="literallyanyletters"
conda activate scanpy-env

python3 xgboost_annotation.py --base_dir "$BASE_DIR" --output_base_dir "$OUTPUT_BASE_DIR"

