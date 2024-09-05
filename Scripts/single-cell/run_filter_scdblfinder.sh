#!/bin/bash

#SBATCH --nodes=1
#SBATCH --output=filter_scdblfinder.out
#SBATCH --error=filter_scdblfinder_error.log
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=6-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --qos=medium


module load anaconda
#export PYTHONNOUSERSITE="literallyanyletters"
conda activate rpy2_env

python3 filter_scdblfinder.py

