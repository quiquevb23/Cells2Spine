'''
    Python file to convert h5ad spatial ref to input required for RCTD
'''

import os
import json
import scanpy as sc
import scipy.sparse as sp
from scipy.io import mmwrite
import pandas as pd

# Load config file
with open("config.json", "r") as f:
    config = json.load(f)

input_dir = config["sc_ref_directory"]
tmp_dir = config["tmp_directory"]

# Ensure tmp directory exists
os.makedirs(tmp_dir, exist_ok=True)

def process_h5ad(adata, label):
    """Extracts counts, cell types, and nUMI from an H5AD file."""
    try:
        output_dir = os.path.join(tmp_dir, label)
        os.makedirs(output_dir, exist_ok=True)  # Create label-specific subdirectory
        
        print(f"Processing {label}...")

        # Convert dense matrix to sparse if necessary
        if not sp.issparse(adata.X):
            print(f"Converting {label} to sparse format")
            adata.X = sp.csr_matrix(adata.X)
        
        # Transpose counts matrix to match RCTD format (genes as rows, cells as columns)
        counts = adata.X.T  

        # Save in Matrix Market (.mtx) format
        counts_path = os.path.join(output_dir, f"counts.mtx")
        mmwrite(counts_path, counts)
        print(f"Saved counts matrix: {counts_path}")

        # Save row (gene) and column (cell) names
        genes_path = os.path.join(output_dir, f"genes.csv")
        cells_path = os.path.join(output_dir, f"cells.csv")
        pd.DataFrame(adata.var_names).to_csv(genes_path, index=False, header=False)
        pd.DataFrame(adata.obs_names).to_csv(cells_path, index=False, header=False)
        print(f"Saved genes: {genes_path}")
        print(f"Saved cells: {cells_path}")

        # Save cell types (if available)
        if 'cell_l4' in adata.obs:
            cell_types_path = os.path.join(output_dir, f"cell_types.csv")
            cell_types = adata.obs[['cell_l4']]
            cell_types.to_csv(cell_types_path, index=True, header=False)  # Keep barcodes as index
            print(f"Saved cell types: {cell_types_path}")
        else:
            print(f"Warning: 'cell_l4' not found in {label}")

        print(f"Finished processing {label}")

    except Exception as e:
        print(f"Error processing {label}: {e}")

if __name__ == "__main__":
    adata_ref = sc.read_h5ad(input_dir)

    # Create subsets based on labels
    subsets = {
        "injured": adata_ref[adata_ref.obs['label'].isin(['7d', '14d'])].copy(),
        "healthy": adata_ref[adata_ref.obs['label'].isin(['uninjured'])].copy()
    }

    for label, adata_subset in subsets.items():
        process_h5ad(adata_subset, label)
    
    print("All H5AD files processed.")