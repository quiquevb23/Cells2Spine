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
from itertools import combinations
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

marker_genes = {
    "Neuron": ["Snap25", "Map2", "Rbfox3", "Syp"],
    "Astrocyte": ["Ntsr2", "Htra1", "Aqp4"],
    "OPC": ["Plp1", "Mobp", "Mag", "Mog"],
    "ODC": ["Gpr17", "Pdgfra", "Sox10"],
    "Microglia": ["Ctss", "Cx3cr1", "Aif1", "Ly86"],
    "Endothelial": ["Cldn5", "Flt1", "Tek", "Cd34","Pecam1", "Prom1"],
    "Pericyte": ["Pdgfrb", "Vtn", "Myl9"],
    "Ependyma": ["Foxj1", "Sox2", "Rsph1", "Ak7"],
    "Stromal": ["Dcn", "Apod", "Gsn", "Col1a1", "Col3a1"],
    "Erythrocyte": ["Hbb-bt", "Hba-a1", "Hba-a2"],
    "Leukocyte": ["Ms4a4b", "Ltb", "Ctsw", "Cd3e"], 
    "Neutrophil": ["S100a8", "S100a9", "Trem1"],
}


# Loop over all subdirectories

adata = sc.read(os.path.join(joined_base_dir, "adata_integrated.h5ad"))

sc.pl.umap(
    adata,
    color=["scDblFinder_score", "scDblFinder_class"],
    wspace=0.5
)
plt.savefig(os.path.join(joined_output_base_dir, 'UMAP_doublet.png'), bbox_inches='tight')
plt.close()

sc.pl.umap(
    adata,
    color=["total_counts", "n_genes_by_counts","pct_counts_mt", "pct_counts_ribo"],
    wspace=0.5
)
plt.savefig(os.path.join(joined_output_base_dir, 'UMAP_QC.png'), bbox_inches='tight')
plt.close()

#Clustering
sc.tl.leiden(adata, flavor="leidenalg", n_iterations=-1, resolution=0.2, key_added="leiden_0.2")
sc.tl.leiden(adata, flavor="leidenalg", n_iterations=-1, resolution=0.4, key_added="leiden_0.4")
sc.tl.leiden(adata, flavor="leidenalg", n_iterations=-1, resolution=0.8, key_added="leiden_0.8")
sc.tl.leiden(adata, flavor="leidenalg", n_iterations=-1, resolution=1.2, key_added="leiden_1.2")

sc.pl.umap(
    adata,
    color=[
        "leiden_0.2",
        "leiden_0.4"
    ],
    wspace=0.5
)
plt.savefig(os.path.join(joined_output_base_dir, 'UMAP_clustering_lowres.png'), bbox_inches='tight')
plt.close()

sc.pl.umap(
    adata,
    color=[
	"leiden_0.8",
        "leiden_1.2"
    ],
    wspace=0.5
)
plt.savefig(os.path.join(joined_output_base_dir, 'UMAP_clustering_hires.png'), bbox_inches='tight')
plt.close()

with open(os.path.join(joined_output_base_dir, 'clustering_info.txt'), 'w') as f:
    f.write(f"Size of adata.raw for joined sample: {adata.raw.shape}\n")
    
    for resolution in ["0.2", "0.4", "0.8", "1.2"]:
        cluster_key = f"leiden_{resolution}"

        nre_cells = adata.obs[cluster_key].value_counts()

        total_cells = len(adata.obs)
        pct_cells = (nre_cells / total_cells) * 100

        sample_1_counts = adata.obs[adata.obs["sample"] == "Sample_1"].groupby(cluster_key).size()
        sample_3_counts = adata.obs[adata.obs["sample"] == "Sample_3"].groupby(cluster_key).size()

        # Fill missing cluster entries with 0
        sample_1_counts = sample_1_counts.reindex(nre_cells.index, fill_value=0)
        sample_3_counts = sample_3_counts.reindex(nre_cells.index, fill_value=0)

        # Calculate percentages for each cluster
        pct_sample_1 = (sample_1_counts / nre_cells) * 100
        pct_sample_3 = (sample_3_counts / nre_cells) * 100

        # Format the percentages as strings and concatenate them
        pct_samples = pct_sample_1.map(lambda x: f"{x:.1f}%") + "/" + pct_sample_3.map(lambda x: f"{x:.1f}%")

        f.write(f"\nSummary statistics for each cluster (resolution {resolution}):\n")
        # Grouping by clusters and computing statistics for total counts and number of genes

        cluster_stats = adata.obs.groupby(cluster_key).agg(
            avg_total_counts=pd.NamedAgg(column="total_counts", aggfunc="mean"),
            avg_genes=pd.NamedAgg(column="n_genes_by_counts", aggfunc="mean"),
            avg_pct_counts_mt=pd.NamedAgg(column="pct_counts_mt", aggfunc="mean"),
            avg_pct_counts_ribo=pd.NamedAgg(column="pct_counts_ribo", aggfunc="mean"),
            doublet_counts=pd.NamedAgg(column="scDblFinder_class", aggfunc=lambda x: (x != 0).sum())
        )

        cluster_stats["nre_cells"] = nre_cells
        cluster_stats["pct_cells"] = pct_cells
        cluster_stats["pct_samples"] = pct_samples

        #Reorder columns
        cluster_stats = cluster_stats[['nre_cells', 'pct_cells', 'pct_samples', 'doublet_counts', 'avg_total_counts', 'avg_genes', 'avg_pct_counts_mt', 'avg_pct_counts_ribo']]

        # Format the DataFrame to ensure it fits well in a text file
        # Optionally, you can set display options for better readability
        pd.set_option('display.max_colwidth', None)  # Ensure no truncation of column widths
        formatted_cluster_stats = cluster_stats.to_string()

        # Writing the cluster statistics directly, with the new column included
        f.write(f"{formatted_cluster_stats}\n")

updated_marker_genes = {}

# Loop through each cell type in the marker_genes dictionary
for cell_type, genes in marker_genes.items():
    # Filter the genes that are present in adata.var_names
    present_genes = [gene for gene in genes if gene in adata.var_names]
    
    # If any genes are present, update the dictionary
    if present_genes:
        updated_marker_genes[cell_type] = present_genes

sc.pl.dotplot(
    adata,
    groupby="leiden_0.2",
    var_names=updated_marker_genes,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    show=False
)
plt.savefig(os.path.join(joined_output_base_dir, 'dotplot_marker_genes_0.2.png'), bbox_inches='tight')
plt.close()

sc.pl.dotplot(
    adata,
    groupby="leiden_0.4",
    var_names=updated_marker_genes,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    show=False
)
plt.savefig(os.path.join(joined_output_base_dir, 'dotplot_marker_genes_0.4.png'), bbox_inches='tight')
plt.close()

sc.pl.dotplot(
    adata,
    groupby="leiden_0.8",
    var_names=updated_marker_genes,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    show=False
)
plt.savefig(os.path.join(joined_output_base_dir, 'dotplot_marker_genes_0.8.png'), bbox_inches='tight')
plt.close()

sc.tl.dendrogram(adata, groupby="leiden_0.2")
sc.tl.rank_genes_groups(adata, groupby="leiden_0.2", method="wilcoxon", use_raw=False)
sc.pl.rank_genes_groups_dotplot(
    adata, groupby="leiden_0.2", standard_scale="var", n_genes=5
)
plt.savefig(os.path.join(joined_output_base_dir, 'dotplot_DEGs_0.2.png'), bbox_inches='tight')
plt.close()

sc.tl.dendrogram(adata, groupby="leiden_0.4")
sc.tl.rank_genes_groups(adata, groupby="leiden_0.4", method="wilcoxon", use_raw=False)
sc.pl.rank_genes_groups_dotplot(
    adata, groupby="leiden_0.4", standard_scale="var", n_genes=5
)
plt.savefig(os.path.join(joined_output_base_dir, 'dotplot_DEGs_0.4.png'), bbox_inches='tight')
plt.close()

sc.tl.dendrogram(adata, groupby="leiden_0.8")
sc.tl.rank_genes_groups(adata, groupby="leiden_0.8", method="wilcoxon", use_raw=False)
sc.pl.rank_genes_groups_dotplot(
    adata, groupby="leiden_0.8", standard_scale="var", n_genes=5
)
plt.savefig(os.path.join(joined_output_base_dir, 'dotplot_DEGs_0.8.png'), bbox_inches='tight')
plt.close()

#If we find as DEGs genes such as Mito, Hb or Ribo, we may need to regress them out
adata.write(os.path.join(joined_base_dir, "joined_DEGs.h5ad"))

'''
# Compute the average expression of marker genes per cluster
for cell_type, genes in marker_genes.items():
    adata.raw.to_adata()
    adata.obs[cell_type] = adata.raw[:, genes].X.mean(axis=1)
    
# Determine dominant cell type per cluster
cluster_marker_expr = adata.groupby('leiden', as_index=False).mean()
dominant_cell_types = cluster_marker_expr.idxmax(axis=1)
adata.obs['cell_type'] = adata.obs['leiden'].map(dominant_cell_types)

cell_types = adata.obs['cell_type'].unique()

# DEA results storage
dea_results = {}


for cell_type in cell_types:
    # Subset data for the specific cell type
    subset = adata[adata.obs['cell_type'] == cell_type]
    
    # Ensure 'sample' column is present for comparison
    if 'sample' not in subset.obs.columns:
        raise ValueError("The 'sample' column is missing from adata.obs.")
    
    # Perform DEA between samples
    for sample1, sample2 in combinations(subset.obs['sample'].unique(), 2):
        subset_sample1 = subset[subset.obs['sample'] == sample1]
        subset_sample2 = subset[subset.obs['sample'] == sample2]
         
        # Perform differential expression analysis
        sc.tl.rank_genes_groups(subset_sample1, groupby='leiden', reference=sample2)
        dea_results[(cell_type, sample1, sample2)] = subset_sample1.uns['rank_genes_groups']
                        
        # Save DEA results
        results_df = pd.DataFrame({
            'gene': dea_results[(cell_type, sample1, sample2)]['names'].flatten(),
            'logfoldchanges': dea_results[(cell_type, sample1, sample2)]['logfoldchanges'].flatten(),
            'pvals': dea_results[(cell_type, sample1, sample2)]['pvals'].flatten(),
            'pvals_adj': dea_results[(cell_type, sample1, sample2)]['pvals_adj'].flatten()
        })
                        
        results_filename = f'dea_{cell_type}_{sample1}_vs_{sample2}.csv'
        results_filepath = os.path.join(output_dir, results_filename)
        results_df.to_csv(results_filepath, index=False)
'''
