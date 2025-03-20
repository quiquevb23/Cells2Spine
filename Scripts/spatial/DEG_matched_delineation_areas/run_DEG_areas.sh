#!/bin/bash

#SBATCH --nodes=1
#SBATCH --output=output_DEG.out
#SBATCH --error=error_DEG.log
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --qos=short


source /home/quiquevb/Cells2Spine/Cells2Spine/Scripts/config_spatial.sh

module load anaconda
#export PYTHONNOUSERSITE="literallyanyletters"
conda activate rpy2_env

python3 DEG_areas.py --base_dir "$BASE_DIR" --output_base_dir "$OUTPUT_BASE_DIR"
