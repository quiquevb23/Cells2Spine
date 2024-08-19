#!/bin/bash

#SBATCH --nodes=1
#SBATCH --output=matrix_output.out
#SBATCH --error=matrix_error.log
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=6-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --qos=medium

export PATH=/storage/gge/Quique/CellRanger/cellranger-8.0.1:$PATH

./count.sh

