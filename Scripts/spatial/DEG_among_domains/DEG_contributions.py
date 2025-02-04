"""
	Script to infer cell type contributions for the DEGs from each cell type:
		-We have a list of DEGs for dorsal region
		-Composition of cell types for that region
		-Single-cell expression data of those cell types

"""


import os
import pandas as pd
from sklearn import metrics
import multiprocessing as mp
import matplotlib.pyplot as plt
import argparse
import scanpy as sc
import scipy.sparse
import numpy as np
import seaborn as sns
import cell2location

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

parent_dir = os.path.dirname(base_dir)
parent_output_dir = os.path.dirname(output_base_dir)

joined_base_dir = os.path.join(parent_dir, "DEG_contribution")
joined_output_base_dir = os.path.join(parent_output_dir, "DEG_contribution")

os.makedirs(joined_base_dir, exist_ok=True)
os.makedirs(joined_output_base_dir, exist_ok=True)


#SC ref
ref_signatures = "/storage/gge/Quique/Cells2SpineData/Pilot/spatial/matrices/Deconvolution/reference_signatures"

#We have here the list for DEGs for dorsal region of injured vs healthy
deg_file_path = "/home/quiquevb/Cells2Spine/Cells2Spine/Outputs/spatial/DEGs/Top_DEG_genes_injured_vs_healthy.csv"

#Celltype proportions
proportions_file_path = "/home/quiquevb/Cells2Spine/Cells2Spine/Outputs/spatial/DEGs/proportions_dorsal_injured.csv"

#Open single-cell ref T.Paralytica and 
def load_references(ref_signatures):
    # load regression model for sc
    adata_file = f"{ref_signatures}/sc.h5ad"
    adata_ref = sc.read_h5ad(adata_file)
    mod = cell2location.models.RegressionModel.load(f"{ref_signatures}", adata_ref)

    adata_ref = mod.export_posterior(
        adata_ref, use_quantiles=True,
        # choose quantiles
        add_to_varm=["q05","q50", "q95", "q0001"],
        #sample_kwargs={'batch_size': 2500, 'use_gpu': False}
    )

    # export estimated expression in each cluster
    if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
        inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_ref.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
    inf_aver.columns = adata_ref.uns['mod']['factor_names']
    return inf_aver

#Load reference signatures from Cell2Location results
cell_type_weights = load_references(ref_signatures)

#Load DEGs
deg_results = pd.read_csv(deg_file_path)
significant_degs = deg_results[
    (deg_results['p_value'] < 0.05) &
    (abs(deg_results['fold_change']) > 1.5)
]


#unique_genes incorporate later (need to order from more to less)

# Get the weight matrix for DEGs and unique genes
# Create a dataframe for weights from the cell type model
weights = cell_type_weights.loc[cell_type_weights.index.isin(deg_results)]

#Add cell-type proportions of region of interest, in our case, dorsal
proportions = pd.read_csv(proportions_file_path, header=None)
proportions = proportions.T  # Transpose the DataFrame
proportions.columns = proportions.iloc[0]  # Set the first row as column names
proportions = proportions[1:]  # Remove the first row from the DataFrame
proportions = proportions.set_index('CellType').astype(float)

# Initialize an empty dataframe to store contributions
contributions = pd.DataFrame(index=weights.index, columns=weights.columns)

# Calculate contributions for each gene
for gene in weights.index:
    contributions.loc[gene] = weights.loc[gene] * proportions

# Step 5: Summarize contributions
# Calculate total contributions across DEGs and unique genes
total_contributions = contributions.sum(axis=0)

# Optionally, visualize the results
total_contributions.plot(kind='bar', figsize=(10, 6))
plt.title('Cell-Type Contributions to DEGs and Unique Genes')
plt.ylabel('Contribution')
plt.xlabel('Cell Types')
plt.xticks(rotation=45)  # Rotate x-axis labels for better visibility
plt.tight_layout()
plt.savefig(os.path.join(joined_output_base_dir("celltype_contributions.png"), bbox_inches='tight')
