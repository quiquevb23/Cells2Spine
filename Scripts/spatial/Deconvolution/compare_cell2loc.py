'''
File to compare cell2loc scores between the following approaches:
    1. Deconvolution with T.Paralytica as reference built from all cells
    2. Deconvolution with T.Paralytica as reference built from either "healthy" or "injured post-7 and post-14 days", separately to
    deconvolute signal for "healthy" and "injured" samples of ST
'''

import os
import pandas as pd
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

base_dir = os.path.join(args.base_dir, "Deconvolution/comparison")
output_dir = os.path.join(args.output_base_dir, "Deconvolution/comparison")

os.makedirs(base_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

# Define paths for files for adata ST data deconvoluted either with separated sc references or single sc

cell2loc_files_separate = '/storage/gge/Quique/Cells2SpineData/Pilot/spatial/matrices/Deconvolution/cell2location_map'
cell2loc_files_together = os.path.join(cell2loc_files_separate, "Deconvolution_wholedataset")
def import_files(path, way):
#Function to import ST data files that have deconvolution scores
    print(way)
    for file in os.listdir(path):
        if file.endswith("h5ad"):
            adata = sc.read(os.path.join(path, file))
            if way == "separate":
                adatas_separate.append(adata)
            elif way == "together":
                adatas_together.append(adata)

def celltype_scores(adatas, way, output_dir):
    print("he")
    num_adatas = len(adatas)
    fig, axes = plt.subplots(num_adatas, 1, figsize=(10, 6 * num_adatas), squeeze=False)

    for i, adata in enumerate(adatas):
        adata.obs[adata.uns['mod']['factor_names']] = adata.obsm['q05_cell_abundance_w_sf']
        name = adata.obs['sample'].unique()[0] # Get name of adata
        
        total_proportions_per_spot = adata.obs[adata.uns['mod']['factor_names']].sum(axis=1)       
        adata.obs['total_proportions'] = total_proportions_per_spot
        
        sc.pl.spatial(adata, color='total_proportions', 
                      title=f"Spatial Plot of Total Cell Type Proportions - {name}", 
                      cmap='viridis')
        plt.savefig(os.path.join(output_dir, f"celltype_per_spot_{name}_{way}.png"))
        plt.close()

        average_proportions = total_proportions_per_spot.mean()
        std_proportions = total_proportions_per_spot.std()

        ax = axes[i][0]
        ax.bar([0], [average_proportions], yerr=[std_proportions], color='skyblue', edgecolor='black', capsize=5)
            
        ax.set_title(f"Total Average Cell Type Proportions: {way} - Dataset {adata.obs['sample'].unique()[0]}", fontsize=16)
        ax.set_ylabel("Proportion", fontsize=14)
        ax.set_xticks([0])
        ax.set_xticklabels(["Total Proportion"], fontsize=14)
        ax.grid(axis='y', linestyle='--', alpha=0.7)
        if way == "separate":
            # Save cell-type proportions for each spot (row-wise) into a CSV file
            celltype_data = adata.obs[adata.uns['mod']['factor_names']]  # This is the dataframe with proportions for each spot
            # Include the cell types (columns) as the header in the CSV file
            output_csv = os.path.join(output_dir, f"celltype_proportions_{name}_{way}.csv")
            celltype_data.to_csv(output_csv)
            print(f"Cell-type proportions saved to {output_csv}")
            
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"celltype_averages_{way}.png"))
    plt.close()
    print("Figure saved to")

# Define the lists to store adatas

adatas_separate = []
adatas_together = []

# Import files with deconvolution results to compute averages for scores

import_files(cell2loc_files_separate, "separate")
import_files(cell2loc_files_together, "together")

celltype_scores(adatas_separate, "separate", output_dir)
celltype_scores(adatas_together, "together", output_dir)