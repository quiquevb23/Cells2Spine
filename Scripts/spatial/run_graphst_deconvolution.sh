#!/bin/bash

#SBATCH --nodes=1
#SBATCH --output=output_graphst_decon.out
#SBATCH --error=error_graphst_decon.log
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --qos=medium


source /home/quiquevb/Cells2Spine/Cells2Spine/Scripts/config_spatial.sh

module load anaconda
#export PYTHONNOUSERSITE="literallyanyletters"
conda activate graphst

python3 graphst_deconvolution.py --base_dir "$BASE_DIR" --output_base_dir "$OUTPUT_BASE_DIR"
