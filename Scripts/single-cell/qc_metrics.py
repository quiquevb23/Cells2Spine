"""
	Script to check QC metrics plots and decide which filters to apply for removing low-quality
	cells
"""
import scanpy as sc
import anndata as ad
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import median_abs_deviation

base_dir = '/storage/gge/Quique/Cells2SpineData/Pilot/single-cell/matrices/'
output_base_dir = '/home/quiquevb/Cells2Spine/Cells2Spine/Outputs/single-cell'

# We will identify cells as outliers if they differ by 5 MADs for the QC in particular

def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name)
    for file in os.listdir(sample_path):
        if file.endswith('filtered_feature_bc_matrix.h5'):
            sample_name = file.replace('_filtered_feature_bc_matrix.h5', '')
            output_dir = os.path.join(output_base_dir, sample_name, 'QC_metrics')
            os.makedirs(output_dir, exist_ok=True)
            adata = sc.read_h5ad(os.path.join(sample_path, file))
            adata.var_names_make_unique()
            adata.layers['counts'] = adata.X #save raw counts
            adata.var['mt'] = adata.var_names.str.startswith('mt-')
            adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
            adata.var["hb"] = adata.var_names.str.startswith(("Hbb", "Hba", "Hbq"))
            adata.raw = adata.copy()
            sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo','hb'], percent_top=[20],
                log1p=True, inplace=True)
            #QC metrics are calculated with unnormalized data and stored, so if we want to use them later they are unnormalized

            sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],
              jitter=0.4, multi_panel=True, show=False)
            plt.savefig(os.path.join(output_dir, 'violin.png'))
            
            sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color='pct_counts_mt', show=False)
            plt.savefig(os.path.join(output_dir, 'scatter.png'))
            
            sns.displot(adata.obs['total_counts'], bins=100, kde=False)
            plt.savefig(os.path.join(output_dir, 'total_counts_displot.png'))
            
            #We will now define if a cell is outlier
            adata.obs["outlier"] = (
                is_outlier(adata, "log1p_total_counts", 5)
                | is_outlier(adata, "log1p_n_genes_by_counts", 5)
                | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
            )
            print(adata.obs.outlier.value_counts())
            #We will filter cells with MADs for pct_counts_mt more 3, and cells with > 8% MT genes
            adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3)
          #      adata.obs["pct_counts_mt"] > 8 #this can be excluded
            print(adata.obs.mt_outlier.value_counts())
            
            adata.write(os.path.join(sample_path, sample_name + "_qc_metrics.h5ad"))

