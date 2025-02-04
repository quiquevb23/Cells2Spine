'''
	Script to normalize individual samples of single-cell and perform PCA and clustering
'''

import numpy as np
import anndata2ri
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import seaborn as sns
import os
import scipy as sp

import rpy2.rinterface_lib.callbacks as rcb
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri

import argparse


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
output_base_dir = args.output_base_dir

ro.pandas2ri.activate()
anndata2ri.activate()

# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name, "outs", "matrices")
    for file in os.listdir(sample_path):
        if file.endswith('filter-low-quality.h5ad'):
            sample_name = file_name

            output_dir = os.path.join(output_base_dir, sample_name)
            os.makedirs(output_dir, exist_ok=True)

            adata = sc.read_h5ad(os.path.join(sample_path, file))
            print(f"Processing sample: {sample_name} with size: {adata.shape}")

            '''
            # Shifted logarithmic normalization
            scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
            adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)
            
            fig, axes = plt.subplots(1, 2, figsize=(10, 5))
            p1 = sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
            axes[0].set_title("Total counts")
            p2 = sns.histplot(adata.layers["log1p_norm"].sum(1), bins=100, kde=False, ax=axes[1])
            axes[1].set_title("Shifted logarithm")
            plt.savefig(os.path.join(output_dir, 'Normalization counts.png'))
            print(adata.shape)
            # HVG selection with scran
            ro.globalenv["adata"] = adata
            '''

            #ro.r('''
            #library(renv)
            #renv::restore()
            #library(scry)
            #library(SingleCellExperiment)
            #sce = devianceFeatureSelection(adata, assay="X") 
            #''')
            
            '''
            binomial_deviance = ro.r("rowData(sce)$binomial_deviance").T
            idx = binomial_deviance.argsort()[-4000:]
            mask = np.zeros(adata.var_names.shape, dtype=bool)
            mask[idx] = True

            adata.var["highly_deviant"] = mask
            adata.var["binomial_deviance"] = binomial_deviance
            '''

            # Scanpy Normalize
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata, n_top_genes=5000, flavor="seurat")
            sc.pl.highly_variable_genes(adata)
            plt.savefig(os.path.join(output_dir, 'HVG_genes.png'))

            if adata.raw is not None:
                print("Adata.raw exists")
#            sc.pp.regress_out(adata, ["pct_counts_mt"])
#            sc.pp.scale(adata)

            adata.write(os.path.join(sample_path, "feature_selection.h5ad"))
