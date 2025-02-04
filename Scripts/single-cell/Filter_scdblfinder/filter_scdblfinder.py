"""
	Script to apply the filtering of low-quality cells and genes, and use scDblFinder in R
"""
import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
from scipy.sparse import issparse
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

# Load R packages and functions
ro.r('''
    library(scDblFinder)
    library(SingleCellExperiment)
''')


# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name, "outs", "matrices")
    for file in os.listdir(sample_path):
        if file.endswith('qc_metrics.h5ad'):
           
            sample_name = file_name
            output_dir = os.path.join(output_base_dir, sample_name)
            os.makedirs(output_dir, exist_ok=True)

            adata = sc.read_h5ad(os.path.join(sample_path, file))
            print(f"Processing sample: {sample_name}")
            sc.pp.filter_genes(adata, min_cells=1) #filter out genes not expressed in any cell
            print(f"Total number of cells: {adata.n_obs}")
            #Filtering cells based on MAD or hard thresholds

            #adata = adata[(~adata.obs.outlier)].copy()
            #We will filter out cells that are outlier for total counts
            adata = adata[~adata.obs['outlier']].copy()
            print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")

            adata = adata[adata.obs['pct_counts_mt'] < 5].copy()
            print(f"Number of cells after filtering those with > 5% pct_mt_counts: {adata.n_obs}")

            adata.raw = adata.copy() #we need to save raw data before normalization but after filtering
            # Convert AnnData object to pandas DataFrame and pass to R
            data_mat = adata.X.T

            if issparse(data_mat):
                data_mat = data_mat.toarray()
                print("Is sparse")
            
            # Count and print the number of NaN values
#            nan_count = np.isnan(data_mat).sum()
#            print(f"Number of NaN values in data_mat: {nan_count}")

#            if nan_count > 0:
#                print("Replacing NaN values with 0")
#                data_mat = np.nan_to_num(data_mat)  # Replace NaN values with 0

#            df_data_mat = pd.DataFrame(data_mat)

            print("DataFrame created from numpy array:")
            
            ro.globalenv['data_mat'] = data_mat
            ro.globalenv['output_dir'] = output_dir
            print("Numpy array successfully converted and transferred to R.")

            # Run scDblFinder in R
            ro.r('''
                print("Reading data as matrix in R")
                set.seed(123)
                print("Converting to SCE and running scDblFinder")
                sce = scDblFinder(
                    SingleCellExperiment(
                        list(counts=data_mat),
                    )
                )

                print("Finish scDblFinder")
                doublet_score <- sce$scDblFinder.score
                doublet_class <- sce$scDblFinder.class
                print("Saving doublet results to output_dir")
                # Save results to R environment
                save(doublet_score, doublet_class, file = file.path(output_dir, "scDblFinder_results.RData"))
                ''')
            print("R_snippet executed")
            # Load results back into Python
            results_path = f"{output_dir}/scDblFinder_results.RData"
            ro.r(f"load('{results_path}')")
            doublet_score = ro.globalenv['doublet_score']
            doublet_class = ro.globalenv['doublet_class']
            print(adata.shape)
            # Convert results back to Python and add to AnnData object
            print("Adding scores and classes to adata")
            adata.obs['scDblFinder_score'] = doublet_score
            adata.obs['scDblFinder_class'] = doublet_class
            print(adata.obs['scDblFinder_class'].value_counts())
            #remove doublets
            adata = adata[adata.obs['scDblFinder_class'] == 1].copy()

            adata.raw = adata.copy() #we need to save raw data before normalization but after filtering
            
            '''
            with open(os.path.join(output_dir, 'doublet_finder_counts.txt'), 'w') as f:
                f.write(f"Doublet finder counts for sample: {sample_name}\n")
                # Write scDblFinder_class counts
                class_counts = adata.obs['scDblFinder_class'].value_counts()
                f.write("scDblFinder_class counts:\n")
                f.write(class_counts.to_string() + "\n\n")
                # Write summary statistics for scDblFinder_score
                score_series = adata.obs['scDblFinder_score']
    
                # Get summary statistics using describe()
                summary = score_series.describe(percentiles=[.25, .5, .75])
    
                # Write the summary statistics to file
                f.write("scDblFinder_score summary statistics:\n")
                f.write(f"Count: {summary['count']}\n")
                f.write(f"Mean: {summary['mean']}\n")
                f.write(f"Std: {summary['std']}\n")
                f.write(f"Min: {summary['min']}\n")
                f.write(f"25th Percentile: {summary['25%']}\n")
                f.write(f"50th Percentile (Median): {summary['50%']}\n")
                f.write(f"75th Percentile: {summary['75%']}\n")
                f.write(f"Max: {summary['max']}\n")
            '''

            # Save the updated AnnData object
            adata.write(os.path.join(sample_path, "filter-low-quality.h5ad"))


            '''
            sc.pp.normalize_total(adata, target_sum=10000)
            sc.pp.log1p(adata)

#            sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=5000)

#            sc.pp.regress_out(adata, ["pct_counts_mt"])
            sc.pp.scale(adata)

            sc.tl.pca(adata, svd_solver="arpack")
            sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
            plt.savefig(os.path.join(output_dir, 'pca_variance.png'))
            adata.write(os.path.join(sample_path, sample_name + "_normalized.h5ad"))
            '''
