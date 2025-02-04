#!/bin/bash

#SBATCH --nodes=1
#SBATCH --output=output_cell2loc.out
#SBATCH --error=error_cell2loc.log
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=15
#SBATCH --time=02-00:00:0
#SBATCH --mem-per-cpu=30G
#SBATCH --qos=medium


source /home/quiquevb/Cells2Spine/Cells2Spine/Scripts/config_spatial.sh

module load anaconda
#export PYTHONNOUSERSITE="literallyanyletters"
conda activate rpy2_env

python3 rctd.py --base_dir "$BASE_DIR" --output_base_dir "$OUTPUT_BASE_DIR"
