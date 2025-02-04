#!/bin/bash

#SBATCH --nodes=1
#SBATCH --output=sctransform.out
#SBATCH --error=sctransform_err.log
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --qos=medium

source /home/quiquevb/Cells2Spine/Cells2Spine/Scripts/config_spatial.sh

module load anaconda
#export PYTHONNOUSERSITE="literallyanyletters"
conda activate rpy2_env

python3 sctransform.py --base_dir "$BASE_DIR" --output_base_dir "$OUTPUT_BASE_DIR"

