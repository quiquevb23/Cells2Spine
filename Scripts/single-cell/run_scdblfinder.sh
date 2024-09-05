#!/bin/bash

#SBATCH --nodes=1
#SBATCH --output=scdblfinder.out
#SBATCH --error=scdblfinder_error.log
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=6-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --qos=medium


eval "$(conda shell.bash hook)"
conda activate rpy2_env

Rscript scdblfinder.R
