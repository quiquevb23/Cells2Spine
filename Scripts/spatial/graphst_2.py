'''
	Script for identification of spatial domains on Visium data using GraphST
	IMPORTANT: it requires log-transformation and normalization of raw counts, scaling and
	selection of top 3000HVGs first
'''
import os
import torch
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp
from GraphST import GraphST
import matplotlib.pyplot as plt

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


# Run device, by default, the package is implemented on 'cpu'. We recommend using GPU.
device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')

# the location of R, which is necessary for mclust algorithm. Please replace the path below with local R installation path
os.environ['R_HOME'] = '/home/quiquevb/.conda/envs/graphst/lib/R'

# the number of clusters
n_clusters = 7


# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name, "outs", "matrices")
    for file in os.listdir(sample_path):
        if file.endswith('qc_metrics.h5ad'):
            sample_name = file_name
            # Define the directories to save plots and matrices

            output_dir = os.path.join(output_base_dir, sample_name)
            os.makedirs(output_dir, exist_ok=True)

            adata = sc.read_h5ad(os.path.join(sample_path, file)) #THIS works but maybe is not getting images

            print(f"Processing sample {sample_name} with size: {adata.shape}")
            
            sc.pp.filter_genes(adata, min_cells=1) #filter out genes not expressed in any cell
            #We will not filter spots for now

            adata.raw = adata.copy() #we need to save raw data before normalization
            
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat")
            
            sc.pl.highly_variable_genes(adata)
            plt.savefig(os.path.join(output_dir, 'HVG_genes.png'))
            plt.close()
            
            adata.write(os.path.join(sample_path, "feature_selection.h5ad"))

            #We select top 3000 HVGs for GraphST model construction
            #It makes more sense to scale AFTER selecting HVGs, but in paper they say to do it before
            adata_hvg = adata[:, adata.var['highly_variable']].copy()  # HVG data for further analysis

            #Let's try 3 approaches here:
            	# 1. Run GraphST with HVG selection
                # 2. Run GraphST with all genes
                # 3. Run GraphST with SVG genes (those selected with Moran's)

            sc.pp.scale(adata_hvg)
            
            #Then train model for spatial domain identification
            model = GraphST.GraphST(adata_hvg, device=device)

            # train model
            adata_hvg = model.train()
            # set radius to specify the number of neighbors considered during refinement
            radius = 20

            tool = 'mclust' # mclust, leiden, and louvain

            # clustering
            from GraphST.utils import clustering

            if tool == 'mclust':
               clustering(adata_hvg, n_clusters, radius=radius, method=tool, refinement=False) # For DLPFC dataset, we use optional refinement step.
            elif tool in ['leiden', 'louvain']:
               clustering(adata_hvg, n_clusters, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=False)
            
            # plotting spatial clustering result
            sc.pl.spatial(adata_hvg,
              img_key="hires",
              color=["domain"],
              show=False)
            plt.savefig(os.path.join(output_dir, 'graphst_domains.png'), bbox_inches='tight')
            plt.close()

            adata_hvg.write(os.path.join(sample_path, "graphst_domains.h5ad"))

