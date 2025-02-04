"""
	Script to Align with PASTE pairwise alignment the Pilot slices of Spatial Transcriptomics, also check if it can be used later for graphst vertical integration
"""


import math
import time
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import style
import paste as pst
import os
import squidpy
import anndata as ad

from skimage import io
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


adatas = {}
# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name, "outs")
    for file in os.listdir(sample_path):
        if file.endswith('filtered_feature_bc_matrix.h5'):
            sample_name = file_name
            adata = squidpy.read.visium(path=sample_path, counts_file=file)
            adata.var_names_make_unique()
            sc.pp.filter_genes(adata, min_cells=1)
            adata.obs["sample"] = file_name
            if file_name == "Spatial_1" or file_name == "Spatial_2":
                adata.obs["condition"] = "healthy"
            else:
                adata.obs["condition"] = "injured"
            i = int(file_name.split('_')[-1])  #Extract number of file
            adatas[i] = adata

start = time.time()
pi12 = pst.pairwise_align(adatas[1], adatas[2])
pi23 = pst.pairwise_align(adatas[2], adatas[3])
pi34 = pst.pairwise_align(adatas[3], adatas[4])

print('Runtime: ' + str(time.time() - start))
pis = [pi12, pi23, pi34]
slices = [adatas[1], adatas[2], adatas[3], adatas[4]]
new_slices, angles, translation  = pst.stack_slices_pairwise(slices, pis, output_params=True, matrix=True)

#New slices is an updated anndata object with all variables and arguments (including alignment coordinates)


#The following plot shows alignment
slice_colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3']

plt.figure(figsize=(7,7))
for i in range(len(new_slices)):
    pst.plot_slice(new_slices[i],slice_colors[i])
plt.legend(handles=[mpatches.Patch(color=slice_colors[0], label='1'),
                    mpatches.Patch(color=slice_colors[1], label='2'),
                    mpatches.Patch(color=slice_colors[2], label='3'),
                    mpatches.Patch(color=slice_colors[3], label='4')])

plt.gca().invert_yaxis()
plt.axis('off')
plt.savefig(os.path.join(joined_output_base_dir,"Aligned_slices.png"), bbox_inches='tight')

for i in range(len(new_slices)):
    print(f"Adata alignment coordinates for slice {i}:")
    print(new_slices[i].obsm["spatial"][:5])
    output_adata_path = os.path.join(joined_base_dir, f"aligned_adata_{i + 1}.h5ad")  # Create a unique filename
    
    new_slices[i].write(output_adata_path)  # Save the AnnData object
    print(f"Saved aligned new_slice_{i + 1} to {output_adata_path}")


'''
# Assuming adata1 and adata2 are your two AnnData objects
adata1 = new_slices[0]  # The first AnnData object
adata2 = new_slices[1]  # The second AnnData object

# Add a new column to `obs` indicating the slice number
adata1.obs['slice'] = 'slice_1'
adata2.obs['slice'] = 'slice_2'

# Concatenate the two AnnData objects
combined_adata = ad.concat([adata1, adata2])

combined_adata.write_h5ad(os.path.join(data_path, "combined_slices.h5ad"))

plt.rcParams["figure.figsize"] = (3, 3)
combined_adata.obsm['spatial'][:, 1] = -1*combined_adata.obsm['spatial'][:, 1]
combined_adata.obs['slice'].replace({'S1':'slice_1', 'S3':'slice_2'}, inplace=True)
ax = sc.pl.embedding(combined_adata, basis='spatial',
                color='slice',
                show=False)
ax.set_title('Aligned image')
#ax.axis('off')
'''
