#!/bin/bash

#SBATCH --nodes=1
#SBATCH --output=output_clustering.out
#SBATCH --error=error_clustering.log
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=00-01:00:0
#SBATCH --mem-per-cpu=5G
#SBATCH --qos=short


source /home/quiquevb/Cells2Spine/Cells2Spine/Scripts/config_spatial.sh

module load anaconda
#export PYTHONNOUSERSITE="literallyanyletters"
conda activate scanpy-env

python3 clustering_deconv.py --base_dir "$BASE_DIR" --output_base_dir "$OUTPUT_BASE_DIR"
