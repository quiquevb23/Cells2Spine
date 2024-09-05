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

base_dir = '/storage/gge/Quique/Cells2SpineData/Pilot/single-cell/matrices/'
output_base_dir = '/home/quiquevb/Cells2Spine/Cells2Spine/Outputs/single-cell/joined'

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
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name)
    for file in os.listdir(sample_path):
        if file.endswith('integrated.h5ad'):
            output_dir = os.path.join(output_base_dir, sample_name, 'QC_metrics')
            sample_name = file.replace('_clustering.h5ad', '')
            adata = sc.read_h5ad(os.path.join(sample_path, file))

sc.tl.leiden(adata, resolution=0.5)  # Adjust resolution parameter as needed

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
