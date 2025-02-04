import os
import re
import scanpy as sc
import pandas as pd
import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
import argparse
import scipy.sparse
import seaborn as sns


cell2loc_file_path = "/storage/gge/Quique/Cells2SpineData/Pilot/spatial/matrices/Deconvolution/cell2location_map"


def load_celltype_prop(adata, sample):
    adata.obs[adata.uns['mod']['factor_names']] = adata.obsm['q05_cell_abundance_w_sf']
    #create df
    df = adata.obs[adata.uns['mod']['factor_names']]
    df_with_indices = df.copy()
   
    int_sample = int(sample[-1]) #get sample name
    adata.obs.index = [f"{idx}_{int_sample}" for idx in adata.obs.index]
    df_with_indices['index'] = adata.obs.index
    
    # Reorder columns to place index first
    df_with_indices = df_with_indices[['index'] + list(df.columns)]
    # Define the file path
    file_path = f"celltype_proportions.csv"
    
    # Check if the file exists
    file_exists = os.path.isfile(file_path)
    
    # Write or append DataFrame to CSV file
    df_with_indices.to_csv(file_path, mode='a', index=False, header=not file_exists)


for file_name in os.listdir(cell2loc_file_path):
    if file_name.startswith("sp"):
        sample = file_name[2:-5]
        adata = sc.read_h5ad(os.path.join(cell2loc_file_path, file_name))
        load_celltype_prop(adata, sample)

