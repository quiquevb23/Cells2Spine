import os
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Define directories
base_dir = '/storage/gge/Quique/Cells2SpineData/Pilot/spatial/'
output_base_dir = '/home/quiquevb/Cells2Spine/Cells2Spine/Outputs/spatial'

# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name)
    for file in os.listdir(sample_path):
        if file.endswith('_qc_metrics.h5ad'):
            output_dir = os.path.join(output_base_dir, sample_name, 'QC_metrics')
            sample_name = file.replace('_qc_metrics.h5ad', '')
            adata = sc.read_h5ad(os.path.join(sample_path, file))
            sc.pp.filter_genes(adata, min_cells=1) #filter out genes not expressed in any cell
            print(f"Total number of cells: {adata.n_obs}")

## 	HERE WE SHOULD APPLY SOME FILTERS FOR REMOVING LOW-QUALITY BARCODES
