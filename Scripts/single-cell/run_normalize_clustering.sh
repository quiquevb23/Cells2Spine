#!/bin/bash

#SBATCH --nodes=1
#SBATCH --output=normalize_clustering.out
#SBATCH --error=normalize_clustering_error.log
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=6-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --qos=medium


module load anaconda
#export PYTHONNOUSERSITE="literallyanyletters"
conda activate scanpy-env

python3 normalize_clustering.py

