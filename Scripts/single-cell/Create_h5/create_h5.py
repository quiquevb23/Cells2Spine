import pandas as pd
import scipy.io
import scipy.sparse
import h5py
import os
import scanpy as sc
import argparse

# Set up argument parser
#parser = argparse.ArgumentParser(description='Process directory argument.')
#parser.add_argument('last_folder', type=str, help='The name of the last folder in the path')

# Parse the arguments
#args = parser.parse_args()

# Define the directory containing the files
directory = '/storage/gge/Quique/TabulaeParalytica/single'
#directory = os.path.join(base_directory, args.last_folder)

# File paths
barcodes_path = os.path.join(directory, 'GSE234774_rnaseq_barcodes.txt')
features_path = os.path.join(directory, 'GSE234774_rnaseq_features.txt')
matrix_path = os.path.join(directory, 'GSE234774_rnaseq_filtered_scRNA.mtx')

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
output_path = os.path.join(directory, 'GSE234774.h5')

# Save AnnData object to HDF5 file
adata.write(output_path)

print("Conversion complete. Data saved to 'GSE234774.h5'.")

