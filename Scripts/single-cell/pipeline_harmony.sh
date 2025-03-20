#!/bin/bash

#source /home/quiquevb/Cells2Spine/Cells2Spine/Scripts/config_cellbender.sh

# Get the absolute path of the current directory
BASE_DIR=$(pwd)

# Set the absolute path for step1, step2, step3
STEP1_DIR="$BASE_DIR/Harmony_integration"
STEP2_DIR="$BASE_DIR/Celltypist"
STEP3_DIR="$BASE_DIR/Ingest_annot"
STEP4_DIR="$BASE_DIR/Compare_Annotations"

#job1=$(sbatch --parsable --chdir=$STEP1_DIR run_harmony_integrate.sh)
#job2=$(sbatch --parsable --dependency=afterok:$job1 --chdir=$STEP2_DIR run_celltypist.sh)
job2=$(sbatch --parsable --chdir=$STEP2_DIR run_celltypist_joined.sh)
job3=$(sbatch --parsable --dependency=afterok:$job2 --chdir=$STEP3_DIR run_ingest_joined.sh)
#job3=$(sbatch --parsable --chdir=$STEP3_DIR run_normalize_HVG.sh)
job4=$(sbatch --parsable --dependency=afterok:$job3 --chdir=$STEP4_DIR run_compare_ann_joined.sh)

