#!/bin/bash

#SBATCH --nodes=1
#SBATCH --output=output.out
#SBATCH --error=error.log
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --qos=short

export PATH=/storage/gge/Quique/CellRanger/cellranger-8.0.1:$PATH

./count.sh

