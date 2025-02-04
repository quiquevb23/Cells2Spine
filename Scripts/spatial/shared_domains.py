
'''
        Script for identification of spatial domains on Visium data using GraphST
        IMPORTANT: it requires log-transformation and normalization of raw counts, scaling and
        selection of top 3000HVGs first (this step is optional)

        Additionally performs Harmony integration of spots to get similar clusters among conditions
'''
import os
import torch
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp
from GraphST import GraphST
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

#Create new folders for "joined" datasets by Harmony
parent_dir = os.path.dirname(base_dir)
parent_output_dir = os.path.dirname(output_base_dir)

joined_base_dir = os.path.join(parent_dir, "joined")
joined_output_base_dir = os.path.join(parent_output_dir, "joined")

os.makedirs(joined_base_dir, exist_ok=True)
os.makedirs(joined_output_base_dir, exist_ok=True)


# Run device, by default, the package is implemented on 'cpu'. We recommend using GPU.
device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')

# the location of R, which is necessary for mclust algorithm. Please replace the path below with loc$
os.environ['R_HOME'] = '/home/quiquevb/.conda/envs/graphst/lib/R'

# the number of clusters
n_clusters = 10

datas_healthy = []
datas_injured = []

# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name, "outs", "matrices")
    for file in os.listdir(sample_path):
        if file.endswith('feature_selection.h5ad'):
            sample_name = file_name

            # Define the directories to save plots and matrices
            output_dir = os.path.join(output_base_dir, sample_name)
            os.makedirs(output_dir, exist_ok=True)

            adata = sc.read_h5ad(os.path.join(sample_path, file))
#            adata.X = adata.raw.X
            #Give raw count matrix or normalized to Harmony Integration or give it normalized?
            # adata will contain a layer [sct_normalized] that we can give as input to Harmony
            adata.obs['sample'] = file_name
            if file_name == "Spatial_1" or file_name == "Spatial_2":
                adata.obs["condition"] = "healthy"
                adata.layers['raw'] = adata.raw.X #store raw in layers for HVG selection
                datas_healthy.append(adata)
            else:
                adata.obs["condition"] = "injured"
                datas_injured.append(adata)

adata_healthy = sc.concat(datas_healthy, index_unique = "_")
adata_injured = sc.concat(datas_injured, index_unique = "_")

#sc.pp.filter_cells(adata_combined, min_counts=1)
#sc.pp.filter_genes(adata_combined, min_counts=1)
#print(f"Size of adata_combined after filtering out cells with 0 counts: {adata_combined.shape}")

'''
# After concatenation of samples, we need to preprocess data again, normalize and everything
adata_combined.layers['counts'] = adata_combined.X.copy() #store raw data
sc.pp.normalize_total(adata_combined) # The right way to normalize count data according to benchmark$
sc.pp.log1p(adata_combined)

# Select top 5000 HVG for combined dataset and filter them
sc.pp.highly_variable_genes(adata_combined, n_top_genes=3000, flavor="seurat")
sc.pl.highly_variable_genes(adata_combined)
plt.savefig(os.path.join(joined_output_base_dir, 'HVG_genes_combined.png'), bbox_inches='tight')
if adata_combined.raw is not None:
    print("Adata.raw exists")
else:
    print("Adata.raw not exists")

# We will not filter out HVG, but rather give them to input for PCA
sc.pp.scale(adata_combined) # Scale before PCA: all genes despite only using HVGs for PCA
'''
# Function to process concatenated objects, with and without integration, and give integrated
# to train model with GraphST
def process_combined(adata_combined):
    # Normalize concatenated object 
    # For now not
    condition = adata_combined.obs['condition'].unique()[0]
    print(f"Size of adata_combined_{condition}: {adata_combined.shape}")
    # Calculate HVGs, needed for GraphST
    sc.pp.highly_variable_genes(adata_combined, n_top_genes=5000, flavor='seurat',inplace=True, batch_key='sample')
    #Seurat v3 paper prioritizes genes that are HVG in both batches, being optimal preprocessing for integrating
    sc.pp.scale(adata_combined) # only scale
    sc.pp.pca(adata_combined, svd_solver='arpack')
    sc.pl.pca_scatter(adata_combined, color="sample")
    plt.savefig(os.path.join(joined_output_base_dir, f"PCA_scatter_{condition}.png"), bbox_inches='tight')
    plt.close()

    sc.pl.pca_variance_ratio(adata_combined)
    plt.savefig(os.path.join(joined_output_base_dir, f"PCA_variance_ratio_{condition}.png"), bbox_inches='tight')
    plt.close()

    # We will compute now UMAP on combined dataset without integration
    adata_noint = adata_combined.copy()
    sc.pp.neighbors(adata_noint)
    sc.tl.leiden(adata_noint, key_added='leiden_1')
    sc.tl.umap(adata_noint)
    sc.pl.umap(adata_noint, color=["sample", "leiden_1"], wspace=0.5)
    plt.savefig(os.path.join(joined_output_base_dir, f"UMAP_unintegrated_{condition}.png"), bbox_inches='tight')
    plt.close()
    adata_noint.write(os.path.join(joined_base_dir, f"adata_no_integrated_{condition}.h5ad"))

    # Now run Harmony
    # Since we have normalized and scaled data after concat, we can run now Harmony, with first 50 PCs (as calculated)
    sc.external.pp.harmony_integrate(adata_combined, key='sample')
    print(adata_combined.obsm['X_pca_harmony'].shape)
    sc.pp.neighbors(adata_combined, use_rep="X_pca_harmony")
    sc.tl.leiden(adata_combined, key_added='leiden_1')
    sc.tl.umap(adata_combined)
    sc.pl.umap(adata_combined, color=["sample", "leiden_1"], wspace=0.5)
    plt.savefig(os.path.join(joined_output_base_dir, f"UMAP_integrated_{condition}.png"), bbox_inches='tight')
    plt.close()
    adata_combined.write(os.path.join(joined_base_dir, f"adata_integrated_{condition}.h5ad"))

    #adata_hvg = adata_combined[:, adata_combined.var['highly_variable']].copy()
    #sc.pp.scale(adata_hvg)
    adata_hvg = adata_combined.copy() #do not filter HVGs
    '''
    #Need to extract sparse matrix to calculate HVGs, mandatory for Graphst
    adata_hvg.layers['raw'] = pd.DataFrame.sparse.from_spmatrix(adata_hvg.raw.X)
    adata_hvg.layers['raw'].columns = adata_hvg.var.index
    adata_hvg.layers['raw'].index = adata_hvg.obs.index
    sc.pp.highly_variable_genes(adata_hvg, n_top_genes=3000, layer='raw', flavor='seurat_v3', inplace=True, batch_key='sample')
    '''

    #Then train model for spatial domain identification
    model = GraphST.GraphST(adata_hvg, device=device)

    # train model
    adata_hvg = model.train()
    # set radius to specify the number of neighbors considered during refinement

    radius = 10
    tool = 'mclust' # mclust, leiden, and louvain

    #Now we won't do refinement step
    # clustering model built only on HVG among 4 samples so we get common clusters
    from GraphST.utils import clustering
    if tool == 'mclust':
        clustering(adata_hvg, n_clusters, radius=radius, method=tool, refinement=True) # $
    elif tool in ['leiden', 'louvain']:
        clustering(adata_hvg, n_clusters, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=True)

    adata_combined.obs["domain"] = adata_hvg.obs["domain"]
    adata_combined.write(os.path.join(joined_base_dir, f"graphst_domains_combined_{condition}.h5ad"))

    # Define a consistent list of clusters (domains) and colors across all samples
    all_domains = adata_hvg.obs['domain'].unique()  # Get all unique domains
    colors = sns.color_palette("tab10", len(all_domains))  # Generate a color palette

    #Then we plot the clustes spatially on each slice
    for file_name in os.listdir(base_dir):
        sample_path = os.path.join(base_dir, file_name, "outs", "matrices")
        for file in os.listdir(sample_path):
            if file.endswith('feature_selection.h5ad'):
                sample_name = file_name
    
                # Define the directories to save plots and matrices
                output_dir = os.path.join(output_base_dir, sample_name)
                os.makedirs(output_dir, exist_ok=True)

                adata = sc.read_h5ad(os.path.join(sample_path, file))
                print(condition)
                print(adata.obs['condition'].unique()[0])
                if condition in adata.obs['condition'].unique():
                    adata.obs['sample'] = file_name
                    filtered_adata_hvg = adata_hvg[adata_hvg.obs['sample'] == sample_name]
                    adata.obs['domain'] = filtered_adata_hvg.obs['domain'].values
                    adata.obs['domain'] = pd.Categorical(adata.obs['domain'], categories=all_domains)
                    adata.obs['leiden_1'] = filtered_adata_hvg.obs['leiden_1'].values
                    # Set the color palette for the 'domain' column
                    adata.uns['domain_colors'] = colors
                    sc.pl.spatial(adata,
                        img_key="hires",
                        color=["domain"],
                        palette=colors,
                        wspace=0.5,
                        show=False)
                    plt.savefig(os.path.join(joined_output_base_dir, f'graphst_domains_{sample_name}.png'),  bbox_inches='tight')
                    plt.close()
                    sc.pl.spatial(adata,
                        img_key="hires",
                        color=["leiden_1"],
                        wspace=0.5,
                        show=False)
                    plt.savefig(os.path.join(joined_output_base_dir, f'leiden_clusters_{sample_name}.png'),  bbox_inches='tight')
                    plt.close()


                    adata.write(os.path.join(sample_path, f"adata_common_domains_{sample_name}.h5ad"))

    sc.pl.umap(adata_combined, color=["sample", "condition", "domain", "leiden_1"], wspace=0.5)
    plt.savefig(os.path.join(joined_output_base_dir, f"domains_umap_{condition}.png"), bbox_inches='tight')
    plt.close()

# Process concatenated samples per condition
process_combined(adata_healthy)
process_combined(adata_injured)
