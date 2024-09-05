'''
	Script to normalize individual samples of single-cell and perform PCA and clustering
'''

import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import os
import scipy as sp

base_dir = '/storage/gge/Quique/Cells2SpineData/Pilot/single-cell/matrices/'
output_base_dir = '/home/quiquevb/Cells2Spine/Cells2Spine/Outputs/single-cell'


# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name)
    for file in os.listdir(sample_path):
        if file.endswith('_filter-low-quality.h5ad'):
            output_dir = os.path.join(output_base_dir, sample_name, 'QC_metrics')
            sample_name = file.replace('_filter-low-quality.h5ad', '')
            adata = sc.read_h5ad(os.path.join(sample_path, file))

            sc.pp.normalize_total(adata, target_sum=10000)
            sc.pp.log1p(adata)

            sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=5000)

#            sc.pp.regress_out(adata, ["pct_counts_mt"])
            sc.pp.scale(adata)
#	Here we store adata files to provide them to Harmony for integration

            adata.write(os.path.join(sample_path, sample_name + "_normalized.h5ad"))

            sc.pp.pca(adata)
            sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
            plt.savefig(os.path.join(output_dir, 'pca_variance.png'))

            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
            sc.tl.leiden(adata, resolution=0.05, key_added="leiden_0.05") #coarse cell types
            sc.tl.leiden(adata, resolution=0.2, key_added="leiden_0.2")
            sc.tl.leiden(adata, resolution=0.5, key_added="leiden_0.5")
            sc.pl.umap(
                adata,
                color=[
                    "leiden_0.05",
                    "leiden_0.2",
                    "leiden_0.5",
                    "scDblFinder_score",
                    "scDblFinder_class",
                    "pct_counts_mt",
                    "total_counts",
                    "pct_counts_hb",
                    "pct_counts_ribo"
                ],
                # Setting a smaller point size to get prevent overlap
                size=2,
            )
            plt.savefig(os.path.join(output_dir, 'UMAP_clustering.png'))
            adata.write(os.path.join(sample_path, sample_name + "_clustering.h5ad"))

