import os
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Define directories
base_dir = '/storage/gge/Quique/Cells2SpineData/Pilot/spatial/'
output_base_dir = '/home/quiquevb/Cells2Spine/Cells2Spine/Outputs/spatial'

# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name)
        for file in os.listdir(sample_path):
            if file.endswith('h5ad'):
                sample_name = file.replace('h5ad', '')
                counts_file = os.path.join(sample_path, file)

        # Create output directories
        output_dir = os.path.join(output_base_dir, sample_name, 'Clustering')
        os.makedirs(output_dir, exist_ok=True)

        # Load the data
        adata = sc.read(counts_file)
        sc.pp.normalize_total(adata, inplace=True)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
        sc.pp.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        sc.tl.leiden(
            adata, key_added="clusters", flavor="igraph", directed=False, n_iterations=2
        )
        plt.rcParams["figure.figsize"] = (4, 4)
        sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4)
        sc.pl.spatial(adata, img_key="hires", color="clusters", size=1.5)
        sc.tl.rank_genes_groups(adata, "clusters", method="t-test")
        sc.pl.rank_genes_groups_heatmap(adata, groups="5", n_genes=10, groupby="clusters")
