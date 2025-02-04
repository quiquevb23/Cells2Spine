#!/bin/bash

#SBATCH --nodes=1
#SBATCH --output=qc_output.out
#SBATCH --error=qc_error.log
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=6-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --qos=medium

source /home/quiquevb/Cells2Spine/Cells2Spine/Scripts/config_spatial.sh

module load anaconda
#export PYTHONNOUSERSITE="literallyanyletters"
conda activate squidpy-env

python3 qc_metrics.py --base_dir "$BASE_DIR" --output_base_dir "$OUTPUT_BASE_DIR"

