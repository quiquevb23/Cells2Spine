#!/bin/bash

#SBATCH --nodes=1
#SBATCH --output=cellbender_3_output.out
#SBATCH --error=cellbender_3_error.log
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=15-00:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --qos=long


module load anaconda
#export PYTHONNOUSERSITE="literallyanyletters"
conda activate cellbender

./cellbender_3.sh

