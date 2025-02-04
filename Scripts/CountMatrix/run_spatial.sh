#!/bin/bash

#SBATCH --nodes=1
#SBATCH --output=matrix_spatial_output.out
#SBATCH --error=matrix_spatial_error.log
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=6-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --qos=medium

export PATH=/storage/gge/Quique/Ranger/SpaceRanger/spaceranger-3.0.1:$PATH

./count_spatial.sh

