'''
	Clustering based on cell2loc results
'''


import os
import pandas as pd
from sklearn import metrics
import multiprocessing as mp
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import scanpy as sc
import scipy.sparse
import numpy as np
import seaborn as sns

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

joined_base_dir = os.path.join(base_dir, "Clustering_Deconv")
joined_output_base_dir = os.path.join(output_base_dir, "Clustering_Deconv")

os.makedirs(joined_base_dir, exist_ok=True)
os.makedirs(joined_output_base_dir, exist_ok=True)


cell2loc_files = '/storage/gge/Quique/Cells2SpineData/Pilot/spatial/matrices/Deconvolution/cell2location_map'

def clustering(adata, name):
    adata.obs[adata.uns['mod']['factor_names']] = adata.obsm['q05_cell_abundance_w_sf']
    proportions = adata.obs[adata.uns['mod']['factor_names']]
    adata_prop = sc.AnnData(X=proportions)
    # Calculate the neighborhood graph
    sc.pp.neighbors(adata_prop, n_neighbors=10, use_rep='X')  # Adjust n_neighbors as needed

    # Perform Leiden clustering
    sc.tl.leiden(adata_prop, resolution=0.4, key_added="ct_leiden_0.4")  # Adjust resolution parameter as needed
    sc.tl.leiden(adata_prop, resolution=0.8, key_added="ct_leiden_0.8")  # Adjust resolution parameter as needed

    # Add cluster labels to the original adata_vis object
    adata.obs['ct_leiden_0.4'] = adata_prop.obs['ct_leiden_0.4']
    adata.obs['ct_leiden_0.8'] = adata_prop.obs['ct_leiden_0.8']

    # Optionally visualize the clustering results
    # Perform UMAP or t-SNE for visualization
    sc.pp.neighbors(adata, n_neighbors=10)
    sc.tl.umap(adata)

    # Plot UMAP with clusters
    sc.pl.umap(adata, color='ct_leiden_0.4', title='Leiden Celltype Proportions 0.4')
    plt.savefig(os.path.join(joined_output_base_dir, f"UMAP_celltype_prop_0.4_{name}.png"), dpi=300, bbox_inches='tight')
    plt.close()
    sc.pl.umap(adata, color='ct_leiden_0.8', title='Leiden Celltype Proportions 0.8')
    plt.savefig(os.path.join(joined_output_base_dir, f"UMAP_celltype_prop_0.8_{name}.png"), dpi=300, bbox_inches='tight')
    plt.close()
    sc.pl.spatial(adata, color=['ct_leiden_0.4', 'ct_leiden_0.8', 'leiden_0.4'])
    plt.savefig(os.path.join(joined_output_base_dir, f"Spatial_celltype_prop_{name}.png"), dpi=300, bbox_inches='tight')
    plt.close()
    cluster_data = pd.DataFrame({
        'leiden_clusters': adata.obs['ct_leiden_0.8'],
        'leiden_0.4': adata.obs['leiden_0.4']
    })

    # Create a crosstab to compare both classifications
    crosstab = pd.crosstab(cluster_data['leiden_clusters'], cluster_data['leiden_0.4'])
    crosstab_col_norm = crosstab.div(crosstab.sum(axis=0), axis=1) * 100

    # Optional: Visualize as a heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(crosstab_col_norm, annot=True, cmap='viridis', fmt='.2f')
    plt.title('Column-normalized correspondence between leiden_clusters and leiden_0.8')
    plt.ylabel('leiden_clusters')
    plt.xlabel('leiden_0.8')
    plt.savefig(os.path.join(joined_output_base_dir, f"Crosstab_leiden_clusters_celltypeprop_clusters_{name}.png"), dpi=300, bbox_inches='tight')
    plt.close()
    return adata

# Function to generate table of cell-type proportions per cluster and sample
def generate_celltype_composition_table(adata, name):
    """
    Generate a table of cell-type compositions (proportions) for each cluster and sample.
    
    Parameters:
        adata : AnnData
            The AnnData object containing the data.
        output_dir : str
            The directory to save the output table.
    """
    factor_names = adata.uns['mod']['factor_names']
    celltype_proportions = adata.obs[factor_names]

    # Define cluster keys
    cluster_keys = ['ct_leiden_0.4', 'ct_leiden_0.8']

    # Prepare output table for each cluster
    for cluster_key in cluster_keys:
        table_data = []
        
        for cluster_id in sorted(adata.obs[cluster_key].unique()):
            # Get data for current cluster
            cluster_data = celltype_proportions[adata.obs[cluster_key] == cluster_id]
            cluster_avg = cluster_data.mean()
            table_data.append(cluster_avg)
        
        # Create a DataFrame from the table data
        table_df = pd.DataFrame(table_data, index=sorted(adata.obs[cluster_key].unique()), columns=factor_names)
        
        # Save table as CSV
        table_df.to_csv(os.path.join(joined_output_base_dir, f"celltype_composition_{cluster_key}_{name}.csv"))

    print(f"Cell-type composition tables saved to {joined_output_base_dir}")


def plot_celltype_proportions(adata, name):
    """
    Plots the proportions of the top 10 most abundant cell types per cluster for each adata object.
    
    Parameters:
        adata : AnnData object
    """    
    # Plot for each adata
    print(f"Processing adata {name}")
        
    factor_names = adata.uns['mod']['factor_names']
    celltype_proportions = adata.obs[factor_names]
    cluster_keys = ['ct_leiden_0.4', 'ct_leiden_0.8']
            
    for cluster_key in cluster_keys:
        # Identify clusters
        
        cluster_ids = sorted(adata.obs[cluster_key].unique())
        fig, axes = plt.subplots(len(cluster_ids), 1, figsize=(10, 6 * len(cluster_ids)), squeeze=False)
        plot_idx = 0    
        for cluster_id in cluster_ids:
            ax = axes[plot_idx, 0]
            # Calculate averages and deviations for the cluster
            cluster_data = celltype_proportions[adata.obs[cluster_key] == cluster_id]
            cluster_avg = cluster_data.mean()
                
            # Select the top 10 cell types for this cluster
            top10_celltypes = cluster_avg.nlargest(10).index
            cluster_avg = cluster_avg[top10_celltypes]
            cluster_std = cluster_data.std()[top10_celltypes]
                
            # Plotting
            x = np.arange(len(top10_celltypes))  # Top 10 cell types
            ax.bar(x, cluster_avg, yerr=cluster_std, alpha=0.8)                
                
            # Formatting
            ax.set_xticks(x)
            ax.set_xticklabels(top10_celltypes, rotation=45, ha='right')
            ax.set_ylabel('Average Proportion')
            ax.set_title(f'Top 10 Cell Types for Cluster {cluster_id} in adata {name}')
            plot_idx += 1  # Move to the next plot in the list
                
        plt.tight_layout()
        plt.savefig(os.path.join(joined_output_base_dir, f"celltype_values_{cluster_key}_{name}.pdf"), bbox_inches='tight')
        plt.close()

adata_list = []

for file in os.listdir(cell2loc_files):
    if file.startswith('sp'):
        # Extract name and prepare output directory
        name = file[2:-5] # Remove 'sp' prefix and '.h5ad' suffix
        output_dir = os.path.join(output_base_dir, name)
        os.makedirs(output_dir, exist_ok=True)
        
        # Load the AnnData object
        adata = sc.read_h5ad(os.path.join(cell2loc_files, file))        
        
        # Perform clustering
        adata = clustering(adata, name)
        
        # Generate the table of cell-type compositions
        generate_celltype_composition_table(adata, name)
        
        plot_celltype_proportions(adata, name)

        # Save the adata files with cell-type clusters
        adata.write(os.path.join(cell2loc_files, f"{name}_ct_clusters.h5ad"))