#!/bin/bash

#SBATCH --nodes=1
#SBATCH --output=output_rctd.out
#SBATCH --error=error_rctd.log
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=15
#SBATCH --time=02-00:00:0
#SBATCH --mem-per-cpu=30G
#SBATCH --qos=medium

# Load config values
CONFIG_FILE="config.json"

INPUT_DIR=$(jq -r '.input_directory' "$CONFIG_FILE")
FILE_TYPES=$(jq -c '.file_types' "$CONFIG_FILE" | tr -d '[]"')
PYTHON_SCRIPTS=$(jq -r '.python_script' "$CONFIG_FILE" | tr -d '[]"' | tr ',' '\n') # Split into individual scripts
R_SCRIPT=$(jq -r '.r_script' "$CONFIG_FILE")
CONDA_ENV=$(jq -r '.conda_env' "$CONFIG_FILE")
RENV_PROJECT=$(jq -r '.renv_project' "$CONFIG_FILE")


source /home/quiquevb/Cells2Spine/Cells2Spine/Scripts/config_spatial.sh

module load anaconda
#export PYTHONNOUSERSITE="literallyanyletters"
echo "Activating Conda environment: $CONDA_ENV"
conda activate rpy2_env

python3 rctd.py --base_dir "$BASE_DIR" --output_base_dir "$OUTPUT_BASE_DIR"
