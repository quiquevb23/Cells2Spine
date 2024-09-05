import os
import numpy as np
import pandas as pd
from scipy.sparse import issparse
import scanpy as sc

base_dir = '/storage/gge/Quique/Cells2SpineData/Pilot/single-cell/matrices/'
output_base_dir = '/home/quiquevb/Cells2Spine/Cells2Spine/Outputs/single-cell'

for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name)
    for file in os.listdir(sample_path):
        if file.endswith('_qc_metrics.h5ad'):
            sample_name = file.replace('_qc_metrics.h5ad', '')
            output_dir = os.path.join(output_base_dir, sample_name, 'QC_metrics')
            adata = sc.read_h5ad(os.path.join(sample_path, file))
            print(f"Processing sample: {sample_name}")

            # Filter out genes not expressed in any cell
            sc.pp.filter_genes(adata, min_cells=1)
            print(f"Total number of cells: {adata.n_obs}")

            # Filter low-quality cells
            adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
            print(f"Number of cells after filtering of low-quality cells: {adata.n_obs}")

            # Extract the expression matrix (adata.X)
            data_mat = adata.X

            # Check if data_mat is sparse, and if so, convert to dense format
            if issparse(data_mat):
                print("Is sparse")
                data_mat = data_mat.toarray()

            # Count and print the number of NaN values
            nan_count = np.isnan(data_mat).sum()
            print(f"Number of NaN values in data_mat: {nan_count}")

            # Replace NaN values with 0 if any are found
            if nan_count > 0:
                print("Replacing NaN values with 0")
                data_mat = np.nan_to_num(data_mat)

            # Convert numpy array to pandas DataFrame
            df_data_mat = pd.DataFrame(data_mat)

            # Save DataFrame to a CSV file for R to read
            csv_path = os.path.join(sample_path, sample_name + "_data_mat.csv")
            df_data_mat.to_csv(csv_path, index=False)

            # Save the filtered AnnData object
            adata.write(os.path.join(sample_path, sample_name + "_filter-low-quality.h5ad"))

            print(f"Data saved to {csv_path}")
            print(f"AnnData object saved to {os.path.join(sample_path, sample_name + '_filter-low-quality.h5ad')}")

