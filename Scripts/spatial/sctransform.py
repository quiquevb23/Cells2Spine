'''
	Script to apply SCT transform from R to normalize spots counts within each slice

'''

import os
import scanpy as sc
import numpy as np
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
from scipy.sparse import issparse
import matplotlib.pyplot as plt
import seaborn as sns
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

pandas2ri.activate()

# Ensure that renv is activated in the R environment
ro.r('''
    library(renv)
    renv::restore()
''')

# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name, "outs", "matrices")
    for file in os.listdir(sample_path):
        if file.endswith('qc_metrics.h5ad'):
            sample_name = file_name
            # Define the directories to save plots and matrices

            output_dir = os.path.join(output_base_dir, sample_name)
            os.makedirs(output_dir, exist_ok=True)

            adata = sc.read_h5ad(os.path.join(sample_path, file)) #THIS works but maybe is not getti$

            print(f"Processing sample {sample_name} with size: {adata.shape}")

            sc.pp.filter_genes(adata, min_cells=1) #filter out genes not expressed in any cell

            adata.raw = adata.copy() #we need to save raw data before normalization
            
            print(adata.X.shape)
            data_mat = adata.X.T

            if issparse(data_mat):
                data_mat = data_mat.toarray()
                print("Is sparse")

            ro.globalenv['data_mat'] = data_mat
            ro.globalenv['output_dir'] = output_dir

            # Run SCTransform in R
            ro.r('''
                library(Seurat)
                library(sctransform)
                set.seed(123)
                
                # Create Seurat Object from data
                counts_matrix <- as.matrix(data_mat)
                seurat_obj <- CreateSeuratObject(counts = counts_matrix)
                
                # Run SCTransform normalization
                seurat_obj <- SCTransform(seurat_obj, verbose = TRUE, min_cells=1)
                
                # Extract normalized data from SCTransform
                normalized_counts <- as.matrix(GetAssayData(seurat_obj, slot = "data"), assay = "SCT")
                
                # Save the normalized data to output_dir
                save(normalized_counts, file = file.path(output_dir, "sctransform_normalized.RData"))
            ''')
            print("R_snippet executed")
            # Load results back into Python
            # Load results back into Python
            results_path = f"{output_dir}/sctransform_normalized.RData"

            ro.r(f"load('{results_path}')")
            normalized_counts = np.array(ro.globalenv['normalized_counts'])

            normalized_counts = normalized_counts.T

            # Add normalized counts back into the AnnData object
            adata.layers["sctransform"] = normalized_counts
            print("SCTransform normalization added to AnnData")

            # Compare normalize total with sctransform
            scanpy_norm = sc.pp.normalize_total(adata, inplace=False) #maybe this step should not be present
            adata.layers['scanpy_norm'] = sc.pp.log1p(scanpy_norm['X'], copy=True)
            # Extract summary statistics from the layers and add them to adata.obs
            adata.obs['sctransform_counts'] = adata.layers['sctransform'].sum(axis=1)  # Sum across genes for each cell
            adata.obs['scanpy_norm_counts'] = adata.layers['scanpy_norm'].sum(axis=1)  # Sum across genes for each cell
            
            plt.rcParams["figure.figsize"] = (8, 8)
            sc.pl.spatial(adata, img_key="hires", color=["total_counts", "sctransform_counts", "scanpy_norm_counts"], wspace=0.5, show=False)
            plt.savefig(os.path.join(output_dir, 'spatial_counts_normalized.png'), bbox_inches='tight')
            plt.close()

            fig, axes = plt.subplots(1, 3, figsize=(10, 5))
            p1 = sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
            axes[0].set_title("Total counts")

            p2 = sns.histplot(adata.layers["sctransform"].sum(1), bins=100, kde=False, ax=axes[1])
            axes[1].set_title("SCT transform")

            p3 = sns.histplot(adata.layers["scanpy_norm"].sum(1), bins=100, kde=False, ax=axes[2])
            axes[2].set_title("Scanpy normalization")

            plt.savefig(os.path.join(output_dir, 'Normalization_counts.png'))
            plt.close()

            adata.write(os.path.join(sample_path, "sctransformed.h5ad"))
            print(f"Processed and saved sample: {sample_name}")

