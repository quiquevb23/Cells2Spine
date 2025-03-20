#!/bin/bash

#SBATCH --nodes=1
#SBATCH --output=output_rctd.out
#SBATCH --error=error_rctd.log
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --time=02-00:00:0
#SBATCH --mem-per-cpu=30G
#SBATCH --qos=medium

# Load config values
CONFIG_FILE="config.json"

INPUT_DIR=$(jq -r '.input_directory' "$CONFIG_FILE")
FILE_TYPES=$(jq -c '.file_types' "$CONFIG_FILE" | tr -d '[]"')
PYTHON_SCRIPTS=$(jq -r '.python_script' "$CONFIG_FILE" | tr -d '[]"' | tr ',' '\n') # Split into individual scripts
R_SCRIPT=$(jq -r '.r_script' "$CONFIG_FILE")
R_CONDA_ENV=$(jq -r '.r_conda_env' "$CONFIG_FILE")
CONDA_ENV=$(jq -r '.conda_env' "$CONFIG_FILE")
RENV_PROJECT=$(jq -r '.renv_project' "$CONFIG_FILE")

# Check if jq is installed
if ! command -v jq &> /dev/null
then
    echo "Error: jq is not installed. Please install jq and try again."
    exit 1
fi

# Activate Conda environment and run Python script
echo "Activating Conda environment: $CONDA_ENV"

module load anaconda
#export PYTHONNOUSERSITE="literallyanyletters"
echo "Activating Conda environment: $CONDA_ENV"
conda activate "$CONDA_ENV"

# Run the first Python script (conversor.py)
echo "Running Python script 1: $(echo $PYTHON_SCRIPTS | cut -d ' ' -f1)"
python3 "$(echo $PYTHON_SCRIPTS | cut -d ' ' -f1)"

# Run the second Python script (conversor_reference)
echo "Running Python script 2: $(echo $PYTHON_SCRIPTS | cut -d ' ' -f2)"
python3 "$(echo $PYTHON_SCRIPTS | cut -d ' ' -f2)"

conda deactivate

conda activate "$R_CONDA_ENV"
# Load R environment and run R script
echo "Running R script: $R_SCRIPT with renv in $RENV_PROJECT"
cd "$RENV_PROJECT"

# Check if renv is initialized
if [ ! -d "renv" ]; then
    echo "Error: renv environment not found in $RENV_PROJECT. Please initialize renv."
    exit 1
fi

Rscript "$R_SCRIPT"

conda deactivate

echo "All processes completed successfully."


