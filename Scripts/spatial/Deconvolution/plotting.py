'''
	Script for better plotting of celltypes deconvolution
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
import re
import anndata

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

joined_base_dir = os.path.join(parent_dir, "Deconvolution")

joined_output_base_dir = os.path.join(parent_output_dir, "Deconvolution/separate_references")
run_name = f'{joined_base_dir}/cell2location_map'

os.makedirs(joined_base_dir, exist_ok=True)
os.makedirs(joined_output_base_dir, exist_ok=True)

# Create placeholders for 4 files
sorted_sample_numbers = [None] * 4
sorted_file_paths = [None] * 4

# Loop through files in the directory and retain only those ending with .h5ad
for filename in os.listdir(f'{joined_base_dir}/cell2location_map'):
    if filename.endswith('.h5ad'):
        # Extract the sample number from the filename using regular expressions
        sample_number = re.search(r'\d+', filename)
        if sample_number:
            sample_number = int(sample_number.group())  # Extract the first number and convert to int
            filepath = os.path.join(f'{joined_base_dir}/cell2location_map', filename)

            # Assign to the correct position in the lists based on the sample number
            # Assuming the sample number ranges from 1 to 4
            if 1 <= sample_number <= 4:
                sorted_sample_numbers[sample_number - 1] = sample_number
                sorted_file_paths[sample_number - 1] = filepath

# Check if all expected files were found
if None in sorted_sample_numbers or None in sorted_file_paths:
    raise ValueError("One or more sample files are missing. Expected sample numbers 1 through 4.")

# Process the sorted files
adata_list = []
for filepath in sorted_file_paths:
    adata_vis = anndata.read_h5ad(filepath)
    adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
    adata_list.append(adata_vis)

# Calculate the global vmax (maximum) for each factor across all samples for homogeneous cmap
global_vmax = {}
for factor in adata_list[0].uns['mod']['factor_names']:
    # Calculate the maximum value across all samples for the current factor
    max_value = max([adata.obs[factor].max() for adata in adata_list])
    global_vmax[factor] = max_value

# Plotting all samples for each cell type on the same figure with consistent color mapping
for factor in adata_list[0].uns['mod']['factor_names']:
    fig, axes = plt.subplots(1, len(adata_list), figsize=(20, 5))  # Subplots for each sample
    for i, adata_vis in enumerate(adata_list):
        ax = axes[i]
        sc.pl.spatial(adata_vis, 
                      cmap='magma', 
                      color=factor, 
                      ncols=1, 
                      size=1.3, 
                      img_key='hires', 
                      vmin=0, 
                      vmax=global_vmax[factor],  # Use the global max value
                      ax=ax,
                      show=False)
        ax.set_title(f"Sample {sorted_sample_numbers[i]}")  # Use the correctly sorted sample number in the title

    # Save the plot for the current cell type
    plt.savefig(os.path.join(joined_output_base_dir, f"deconvolution_{factor}.png"), dpi=300)
    plt.close()
