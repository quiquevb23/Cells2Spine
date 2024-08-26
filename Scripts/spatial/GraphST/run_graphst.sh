#!/bin/bash

#SBATCH --nodes=1
#SBATCH --output=output_graphst.out
#SBATCH --error=error_graphst.log
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=6-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --qos=medium


module load anaconda
#export PYTHONNOUSERSITE="literallyanyletters"
conda activate graphst

python3 graphst.py

