'''
	Script for dimensionality reduction and clustering of individual samples of single-cell
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

for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name, "outs", "matrices")
    for file in os.listdir(sample_path):
        if file.endswith('feature_selection.h5ad'):
            sample_name = file_name

            output_dir = os.path.join(output_base_dir, sample_name)

            adata = sc.read_h5ad(os.path.join(sample_path, file))

            #Dimenstionality reduction
            #HVG are set as highly deviant (by scry), and selected to compute PCA only on them
#            adata.var["highly_variable"] = adata.var["highly_deviant"]
#            adata.X = adata.layers["log1p_norm"] #we get normalized data
            #Compute PCA using HVG, no need to filter out HVG
            #We need to scale data before running PCA, since Highly expressed genes can mask out PCA resutls

            sc.pp.scale(adata)
            sc.pp.pca(adata, svd_solver="arpack", mask_var="highly_variable")

            sc.pl.pca_scatter(adata, color="total_counts")
            plt.savefig(os.path.join(output_dir, "PCA_scatter.png"))
            plt.close()
            sc.pl.pca_variance_ratio(adata)
            plt.savefig(os.path.join(output_dir, "PCA_variance_ratio.png"))
            plt.close()

            #We store adata here to give to Harmony for integration, with PCs calculated but no neighborhood graph
            adata.write(os.path.join(sample_path, "PCAs.h5ad"))
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)

            #Rename scDblFinder_class as categorical to correctly plot
            adata.obs["scDblFinder_class"] = pd.Categorical(adata.obs["scDblFinder_class"])
            sc.pl.umap(
                adata,
                wspace=0.5,
                color=["scDblFinder_score", "scDblFinder_class"],
            )
            plt.savefig(os.path.join(output_dir, "UMAP_doublet.png"))
            plt.close()

            sc.pl.umap(
                adata,
                wspace=0.5,
                color=["total_counts", "n_genes_by_counts","pct_counts_mt", "pct_counts_ribo"],
            )
            plt.savefig(os.path.join(output_dir, "UMAP_QC.png"))
            plt.close()

            #Clustering
            sc.tl.leiden(adata, flavor="leidenalg", n_iterations=-1, resolution=0.2, key_added="leiden_0.2")
            sc.tl.leiden(adata, flavor="leidenalg", n_iterations=-1, resolution=0.4, key_added="leiden_0.4")
            sc.tl.leiden(adata, flavor="leidenalg", n_iterations=-1, resolution=0.8, key_added="leiden_0.8")
            sc.tl.leiden(adata, flavor="leidenalg", n_iterations=-1, resolution=1.2, key_added="leiden_1.2")

            sc.pl.umap(
                adata,
                wspace=0.5,
                color=[
                    "leiden_0.2",
                    "leiden_0.4",
                ]
            )
            plt.savefig(os.path.join(output_dir, 'UMAP_clustering_lowres.png'))
            plt.close()

            sc.pl.umap(
                adata,
                wspace=0.5,
                color=[
                    "leiden_0.8",
                    "leiden_1.2",
                ]
            )
            plt.savefig(os.path.join(output_dir, 'UMAP_clustering_hires.png'))
            plt.close()

            print(f"Size of adata.raw for sample: {sample_name}: {adata.raw.shape}")
            with open(os.path.join(output_dir, 'clustering_info.txt'), 'w') as f:
                f.write(f"Size of adata.raw for sample {sample_name}: {adata.raw.shape}\n")
    
                for resolution in ["0.2", "0.4", "0.8", "1.2"]:
                    cluster_key = f"leiden_{resolution}"

                    nre_cells = adata.obs[cluster_key].value_counts()
                    total_cells = len(adata.obs)
                    pct_cells = (nre_cells / total_cells) * 100

                    f.write(f"\nSummary statistics for each cluster (resolution {resolution}):\n")
        
                    # Grouping by clusters and computing statistics for total counts and number of genes
                    cluster_stats = adata.obs.groupby(cluster_key).agg(
                        avg_total_counts=pd.NamedAgg(column="total_counts", aggfunc="mean"),
                        avg_genes=pd.NamedAgg(column="n_genes_by_counts", aggfunc="mean"),
                        avg_pct_counts_mt=pd.NamedAgg(column="pct_counts_mt", aggfunc="mean"),
                        avg_pct_counts_ribo=pd.NamedAgg(column="pct_counts_ribo", aggfunc="mean"),
                        doublet_counts=pd.NamedAgg(column="scDblFinder_class", aggfunc=lambda x: (x != 0).sum()),
                    )
                    cluster_stats["nre_cells"] = nre_cells
                    cluster_stats["pct_cells"] = pct_cells
                    #Reorder columns
                    cluster_stats = cluster_stats[['nre_cells', 'pct_cells', 'doublet_counts', 'avg_total_counts', 'avg_genes', 'avg_pct_counts_mt', 'avg_pct_counts_ribo']]
                    # Optionally, you can set display options for better readability
                    pd.set_option('display.max_colwidth', None)  # Ensure no truncation of column widths
                    formatted_cluster_stats = cluster_stats.to_string()
                    f.write(f"{formatted_cluster_stats}\n")

            adata.write(os.path.join(sample_path, "clustering.h5ad"))

