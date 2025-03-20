#!/bin/bash

#SBATCH --nodes=1
#SBATCH --output=output_clusterProf.out
#SBATCH --error=error_clusterProf.log
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --qos=short

# Define config file (or pass as argument)
CONFIG_FILE=${1:-config.json}

# Extract config values
R_SCRIPT=$(jq -r '.r_script' "$CONFIG_FILE")
RENV_PROJECT=$(jq -r '.renv_project' "$CONFIG_FILE")


# Create a new conda environment with R 4.4
ENV_NAME="r_4.4_env"

module load anaconda
# Activate the environment
conda activate "$ENV_NAME"

# Ensure renv project directory is correct
echo "Restoring renv environment in project: $RENV_PROJECT"
Rscript -e "renv::restore(project = '$RENV_PROJECT')"

# Run the R script
echo "Running R script: $R_SCRIPT"
Rscript "$R_SCRIPT"

echo "Functional Enrichment done."

