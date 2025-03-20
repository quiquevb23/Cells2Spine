#!/bin/bash

#source /home/quiquevb/Cells2Spine/Cells2Spine/Scripts/config_cellbender.sh

# Get the absolute path of the current directory
BASE_DIR=$(pwd)

# Set the absolute path for step1, step2, step3
STEP1_DIR="$BASE_DIR/Qc_metrics"
STEP2_DIR="$BASE_DIR/Filter_scdblfinder"
STEP3_DIR="$BASE_DIR/Normalize_HVG"
STEP4_DIR="$BASE_DIR/Clustering"
STEP5_DIR="$BASE_DIR/Celltypist"
STEP6_DIR="$BASE_DIR/Ingest_annot"
STEP7_DIR="$BASE_DIR/Compare_Annotations"

#job1=$(sbatch --parsable --chdir=$STEP1_DIR run_qc_metrics.sh)
#job2=$(sbatch --parsable --dependency=afterok:$job1 --chdir=$STEP2_DIR run_filter_scdblfinder.sh)
job2=$(sbatch --parsable --chdir=$STEP2_DIR run_filter_scdblfinder.sh)
job3=$(sbatch --parsable --dependency=afterok:$job2 --chdir=$STEP3_DIR run_normalize_HVG.sh)
#job3=$(sbatch --parsable --chdir=$STEP3_DIR run_normalize_HVG.sh)
job4=$(sbatch --parsable --dependency=afterok:$job3 --chdir=$STEP4_DIR run_clustering.sh)
job5=$(sbatch --parsable --dependency=afterok:$job4 --chdir=$STEP5_DIR run_celltypist.sh)
#job5=$(sbatch --parsable --chdir=$STEP5_DIR run_celltypist.sh)
job6=$(sbatch --parsable --dependency=afterok:$job5 --chdir=$STEP6_DIR run_ingest.sh)
#job6=$(sbatch --parsable --chdir=$STEP6_DIR run_ingest.sh)
job7=$(sbatch --parsable --dependency=afterok:$job6 --chdir=$STEP7_DIR run_compare_ann.sh)


