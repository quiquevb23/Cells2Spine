#!/bin/bash

#SBATCH --nodes=1
#SBATCH --output=output_graphst_2.out
#SBATCH --error=error_graphst_2.log
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=6-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --qos=medium


source /home/quiquevb/Cells2Spine/Cells2Spine/Scripts/config_spatial_2.sh

module load anaconda
#export PYTHONNOUSERSITE="literallyanyletters"
conda activate graphst

python3 graphst_2.py --base_dir "$BASE_DIR" --output_base_dir "$OUTPUT_BASE_DIR"
