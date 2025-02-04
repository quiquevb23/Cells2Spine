#!/bin/bash

#SBATCH --nodes=1
#SBATCH --output=cellbender_4_output.out
#SBATCH --error=cellbender_4_error.log
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=15-00:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --qos=long


module load anaconda
#export PYTHONNOUSERSITE="literallyanyletters"
conda activate cellbender

./cellbender_4.sh

