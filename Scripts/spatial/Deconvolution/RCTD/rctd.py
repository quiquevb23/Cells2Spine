'''
    Script to perform deconvolution with our own annotated data
'''

import os
import pandas as pd
from sklearn import metrics
import multiprocessing as mp
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import scanpy as sc
import scipy.sparse
import numpy as np
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from scipy.sparse import issparse

def parse_args():
    parser = argparse.ArgumentParser(description="Process directories for single-cell data.")
    
    # Define the arguments for base_dir and output_base_dir
    parser.add_argument('--base_dir', type=str, required=True, 
                        help="Base directory for input data.")
    parser.add_argument('--output_base_dir', type=str, required=True, 
                        help="Base directory for output data.")
    
    # Parse the arguments
    return parser.parse_args()

args = parse_args()

base_dir = args.base_dir
parent_dir = os.path.dirname(base_dir)
model_dir = os.path.join(parent_dir, "celltypist_models")
os.makedirs(model_dir, exist_ok=True)

output_base_dir = args.output_base_dir

# Ensure that renv is activated in the R environment
ro.r('''
    library(renv)
    renv::restore()
''')

# Load R packages and functions
ro.r('''
    library(SingleR)
    library(SummarizedExperiment)
    run_singleR <- function(query_matrix, ref_matrix, ref_labels) {
        # Convert query and reference matrices to SummarizedExperiment
        query_se <- SummarizedExperiment(assays = list(counts = query_matrix))
        ref_se <- SummarizedExperiment(assays = list(counts = ref_matrix))

        # Run SingleR
        results <- SingleR(test = query_se, ref = ref_se, labels = ref_labels)

        # Return annotations
        return(data.frame(Cell = rownames(results), Labels = results$labels))
    }
''')

adata_dict = {}

def deconvolution(adata, sample_name, condition):
    """
        Do the deconvolution with RCTD
    """


# Load references from annotated single-cell data
#adata_ref = 

# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name, "outs", "matrices")
    for file in os.listdir(sample_path):
        if file.startswith('adata_stlearn_common_domains_'):
            sample_name = file_name
            # Define the directories to save plots and matrices
            output_dir = os.path.join(output_base_dir, sample_name)
            os.makedirs(output_dir, exist_ok=True)
            # This has the annotation for clusters but is not normalized
            adata = sc.read_h5ad(os.path.join(sample_path, file))
            condition = adata.obs['condition'].unique()[0]
            deconvolution(adata, sample_name, condition)

