#!/bin/bash

#Script to aggregate matrices from single-cell

#We need to have samples stored in the following directories, with sample name subfolders and 
#fastqfiles subfolder in each

cellranger aggr --id=Sample_1 --csv=./aggregate_csv/sample1_fullseq_aggr.csv
cellranger aggr --id=Sample_3 --csv=./aggregate_csv/sample3_fullseq_aggr.csv

echo "Aggregate completed."



