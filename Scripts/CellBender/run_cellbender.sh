#!/bin/bash

#SBATCH --nodes=1
#SBATCH --output=cellbender_output.out
#SBATCH --error=cellbender_error.log
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=6-00:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --qos=medium


module load anaconda
#export PYTHONNOUSERSITE="literallyanyletters"
conda activate cellbender

./cellbender.sh

