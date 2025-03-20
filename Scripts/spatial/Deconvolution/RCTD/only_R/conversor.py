'''
    Python file to convert h5ad spatial h5ad files to input required for RCTD
'''

import os
import json
import scanpy as sc
import scipy.sparse as sp
import pandas as pd
import multiprocessing


# Load config file
with open("config.json", "r") as f:
    config = json.load(f)

input_dir = config["input_directory"]
tmp_dir = config["tmp_directory"]
num_processes = config["num_processes"]

# Ensure tmp directory exists
os.makedirs(tmp_dir, exist_ok=True)
print(f"Writing output to: {tmp_dir}")

def process_h5ad(file_path):
    print(f"Reading: {file_path}")
    """Extracts counts, cell types, and nUMI from an H5AD file."""
    try:
        # Get the filename from the full path  
        file_name = os.path.basename(file_path)  # Extracts "adata_Spatial_2_manual_delineation.h5ad"
        
        # Extract "Spatial_2" from "adata_Spatial_2_manual_delineation.h5ad"
        parts = file_name.split("_")
        sample_name = "_".join(parts[1:3]) if len(parts) > 2 else file_name  # "Spatial_2"

        print(f"Processing {sample_name}...")
        
        adata = sc.read_h5ad(file_path)
        print(f"Successfully loaded {file_path}")
        # Extract and transpose the counts matrix (sparse format)
        if 'counts' in adata.layers:
            counts = adata.layers['counts'].T  # Transpose: genes (rows) x spots (cols)
        elif sp.issparse(adata.X):  # Check if adata.X is already sparse
            print(f"Using sparse adata.X for {file_name}")
            counts = adata.X.T
        else:
            print(f"Converting dense adata.X to sparse format for {file_name}")
            adata.X = sp.csr_matrix(adata.X)  # Convert to sparse format
            counts = adata.X.T
        
        counts_df = pd.DataFrame.sparse.from_spmatrix(
            counts, index=adata.var_names, columns=adata.obs_names
        )
        try:
            output_path = os.path.join(tmp_dir, f"{sample_name}_counts.csv")
            counts_df.to_csv(output_path)
            print(f"Saved: {output_path}")
        except Exception as e:
            print(f"Error writing {output_path}: {e}")


        # Extract coordinates (array_row and array_col) from .obs
        if 'array_row' in adata.obs and 'array_col' in adata.obs:
            coords = adata.obs[['array_row', 'array_col']]
            coords.to_csv(os.path.join(tmp_dir, f"{sample_name}_coordinates.csv"))
        else:
            print(f"Warning: 'array_row' or 'array_col' not found in {file_name}")

        print(f"Finished processing {sample_name}")
    except Exception as e:
        print(f"Error processing {file_path}: {e}")

if __name__ == "__main__":
    # Find all H5AD files in the input directory
    
    h5ad_files = []
    for file_name in os.listdir(os.path.join(input_dir, "indiv_samples")):
        sample_path = os.path.join(input_dir, "indiv_samples", file_name, "outs", "matrices")
        if os.path.isdir(sample_path):  # Ensure it's a directory
            for file in os.listdir(sample_path):
                if file.endswith("manual_delineation.h5ad"):
                    h5ad_files.append(os.path.join(sample_path, file))
    print(f"Found {len(h5ad_files)} H5AD files: {h5ad_files}")
    for file in h5ad_files:
        process_h5ad(file)
    '''
    # Process files in parallel
    with multiprocessing.Pool(num_processes) as pool:
        pool.map(process_h5ad, h5ad_files)
    '''
    print("All H5AD files processed.")

