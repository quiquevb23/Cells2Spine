'''
	Script for doing GraphST vertical integration based on condition:
	find common domains based on condition
	It needs to be run on the aligned images with PASTE
'''
import os
import torch
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp
from GraphST import GraphST
import matplotlib.pyplot as plt
import anndata as ad
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

#Create new folders for "joined" datasets by Harmony
parent_dir = os.path.dirname(base_dir)
parent_output_dir = os.path.dirname(output_base_dir)

joined_base_dir = os.path.join(parent_dir, "paste_aligned")
joined_output_base_dir = os.path.join(parent_output_dir, "paste_aligned")

os.makedirs(joined_base_dir, exist_ok=True)
os.makedirs(joined_output_base_dir, exist_ok=True)


# Run device, by default, the package is implemented on 'cpu'. We recommend using GPU.
device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')

# the location of R, which is necessary for mclust algorithm. Please replace the path below with loc$
os.environ['R_HOME'] = '/home/quiquevb/.conda/envs/graphst/lib/R'

# the number of clusters
n_clusters = 7

adata_healthy = []
adata_injured = []
# Loop over all subdirectories
for file_name in os.listdir(joined_base_dir):
    if file_name.startswith('aligned'):
        file_path = os.path.join(joined_base_dir, file_name)  # Create the full file path
        adata = sc.read_h5ad(file_path)  # Read the AnnData file
        condition = adata.obs['condition'].unique()[0]
        if condition == "healthy":
            adata_healthy.append(adata)
        else:
            adata_injured.append(adata)

if len(adata_healthy) == 2:
    combined_healthy = ad.concat(adata_healthy, label='sample', join='outer')  # Concatenate healthy samples
    combined_healthy_path = os.path.join(joined_base_dir, "combined_healthy_adata.h5ad")
    combined_healthy.var_names_make_unique()
    combined_healthy.write(combined_healthy_path)  # Save the combined AnnData
    print(f"Saved the combined AnnData object for condition 'healthy' to {combined_healthy_path}")
else:
    print(f"Skipping concatenation for 'healthy': expected 2 samples, found {len(adata_healthy)}.")

if len(adata_injured) == 2:
    combined_injured = ad.concat(adata_injured, label='sample', join='outer')  # Concatenate injured samples
    combined_injured_path = os.path.join(joined_base_dir, "combined_injured_adata.h5ad")
    combined_injured.var_names_make_unique()
    combined_injured.write(combined_injured_path)  # Save the combined AnnData
    print(f"Saved the combined AnnData object for condition 'injured' to {combined_injured_path}")
else:
    print(f"Skipping concatenation for 'injured': expected 2 samples, found {len(adata_injured)}.")

### Get original images
for file_name in os.listdir(base_dir):
    if file_name == "Spatial_1":
        sample_path = os.path.join(base_dir, file_name, "outs", "matrices")
        for file in os.listdir(sample_path):
            if file.endswith('qc_metrics.h5ad'):
                adata = sc.read_h5ad(os.path.join(sample_path, file))
                img_key_healthy = adata.uns["spatial"] #get the hires image coords (probably have to be edited with output of PASTE)
    elif file_name == "Spatial_3":
        sample_path = os.path.join(base_dir, file_name, "outs", "matrices")
        for file in os.listdir(sample_path):
            if file.endswith('qc_metrics.h5ad'):
                adata = sc.read_h5ad(os.path.join(sample_path, file))
                img_key_injured = adata.uns["spatial"] #get the hires image coords (probably have to be $


### RUN GraphST

def plt_aligned_image(adata):
    plt.rcParams["figure.figsize"] = (3, 3)
    condition = adata.obs['condition'].unique()[0]
    adata.obsm['spatial'][:, 1] = -1*adata.obsm['spatial'][:, 1] #This flips the Y-axis
    sc.pl.embedding(adata, basis='spatial',
                    color='sample',
                    show=False)
    plt.savefig(os.path.join(joined_output_base_dir, f'aligned_embeddings_{condition}.png'), bbox_inches='tight')
    plt.close()
    return

def run_graphst(adata):
    condition = adata.obs['condition'].unique()[0]
    # define model
    model = GraphST.GraphST(adata, device=device)
    # run model
    adata = model.train()
    # clustering
    from GraphST.utils import clustering

    tool = 'mclust' # mclust, leiden, and louvain

    # clustering
    from GraphST.utils import clustering

    if tool == 'mclust':
       clustering(adata, n_clusters, method=tool, radius=16, refinement=True) # For DLPFC dataset, we use optional refinement step.
    elif tool in ['leiden', 'louvain']:
       clustering(adata, n_clusters, method=tool, start=0.1, end=2.0, increment=0.01)
    
    ### Plotting UMAP before batch effect correction
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.pca(adata)

    sc.pp.neighbors(adata, use_rep='X_pca', n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color='sample', title='Uncorrected',
                  show=False)
    plt.savefig(os.path.join(joined_output_base_dir, f'UMAP_uncorrected_{condition}.png'), bbox_inches='tight')
    plt.close()

    ### Plotting UMAP after batch effect correction
    sc.pp.neighbors(adata, use_rep='emb_pca', n_neighbors=10)
    sc.tl.umap(adata)
    sc.pl.umap(adata,
           color='sample',
           title='Batch corrected',
           #legend_loc = 'bottom margin',
           show=False)
    plt.savefig(os.path.join(joined_output_base_dir, f'UMAP_batch_corrected_{condition}.png'), bbox_inches='tight')
    plt.close()

    ### Color by predicted domains
    sc.pl.umap(adata, color='domain', title='Colored by clusters', show=False)
    plt.savefig(os.path.join(joined_output_base_dir, f'domains_spatial_{condition}.png'), bbox_inches='tight')
    plt.close()

    plt.tight_layout(w_pad=0.02)
    
    adata.uns['spatial'] = img_key_healthy if condition == "healthy" else img_key_injured
    # plotting spatial clustering result
    sc.pl.spatial(adata,
                  img_key="hires",
                  color=["domain"],
                  show=False)
    plt.savefig(os.path.join(joined_output_base_dir, f'graphst_domains_{condition}.png'), bbox_inches='tight')
    plt.close()
    adata.write(os.path.join(joined_base_dir, f"graphst_domains_{condition}.h5ad"))
    return

plt_aligned_image(combined_healthy)
plt_aligned_image(combined_injured)

run_graphst(combined_healthy)
run_graphst(combined_injured)

