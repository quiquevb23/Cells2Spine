import pandas as pd
import scipy.io
import scipy.sparse
import h5py
import os
import scanpy as sc

# Define the directory containing the files
directory = '/storage/gge/Quique/Cells2SpineData/Pilot/single-cell/matrices/sci_mouse'

# File paths
barcodes_path = os.path.join(directory, 'barcodes.tsv.gz')
features_path = os.path.join(directory, 'features.tsv.gz')
matrix_path = os.path.join(directory, 'matrix.mtx.gz')

# Read barcodes
barcodes_df = pd.read_csv(barcodes_path, sep='\t', header=None)
barcodes = barcodes_df[0].values.astype('U')  # Convert to numpy array and byte strings

# Read features
features_df = pd.read_csv(features_path, sep='\t', header=None)
features = features_df.values.astype('U')  # Convert to numpy array and byte strings

# Read matrix
matrix = scipy.io.mmread(matrix_path)

# Convert matrix to dense format if it is in sparse format
if scipy.sparse.issparse(matrix):
    matrix = matrix.toarray()

matrix = matrix.T

# Create an AnnData object
adata = sc.AnnData(
    X=matrix,
    obs=pd.DataFrame(index=barcodes),
    var=pd.DataFrame(index=features[:, 1])  # Assuming feature names are in the second column
)

# Define the output path
output_path = os.path.join(directory, 'sci_mouse_filtered_feature_bc_matrix.h5')

# Save AnnData object to HDF5 file
adata.write(output_path)

print("Conversion complete. Data saved to 'data.h5'.")

