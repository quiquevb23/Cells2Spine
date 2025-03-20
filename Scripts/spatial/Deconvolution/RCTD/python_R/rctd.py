'''
    Script to perform deconvolution with our own annotated data
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
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from scipy.sparse import issparse

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
parent_dir = os.path.dirname(base_dir)
model_dir = os.path.join(parent_dir, "celltypist_models")
os.makedirs(model_dir, exist_ok=True)

output_base_dir = args.output_base_dir

#Create new folders for "joined" datasets by Harmony
parent_dir = os.path.dirname(base_dir)
parent_output_dir = os.path.dirname(output_base_dir)

#ref_signatures = os.path.join(joined_base_dir, "reference_signatures")
output_dir = os.path.join(parent_output_dir, "Deconvolution", "RCTD")

#os.makedirs(ref_signatures, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

pandas2ri.activate()

# Ensure that renv is activated in the R environment
ro.r('''
    library(renv)
    renv::restore()
''')

def make_ref(adata, directory):
    # Get condition
    condition = 'healthy' if 'uninjured' in adata.obs['label'].values else 'injured'
    
    counts = pd.DataFrame(adata.X.T) # Transpose to have genes as rows and cells as columns
    cell_types = pd.Series(adata.obs['cell_l4'])  # Get layer 4 column exists
    nUMI = pd.Series(counts.sum(axis=0))  # Total UMI counts per cell (sum across column)
    
    # Convert to R objects
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_counts = ro.conversion.py2rpy(pd.DataFrame(counts))  # Convert counts DataFrame to R
        r_cell_types = ro.conversion.py2rpy(pd.Series(cell_types))  # Convert cell types to R
        r_nUMI = ro.conversion.py2rpy(pd.Series(nUMI))  # Convert nUMI to R

    r_script = f"""    
    reference <- Reference({r_counts}, {r_cell_types}, {r_nUMI})
    saveRDS(reference, "{directory}/SCRef_{condition}.rds")
    """
    
    ro.r(r_script)

def deconvolution(adata, sample_name, condition, sample_output_dir, ref_dir):
    """
        Do the deconvolution with RCTD
    """
    ref_file = os.path.join(ref_dir, f"SCRef_{condition}.rds")
    
    # coords need to be rows for each barcode and 2 cols for X and Y (spatial pixels)

    coords = adata.obs[['array_row', 'array_col']].copy()
    coords.columns = ['x', 'y'] # rename for consistency

    spatial_counts = pd.DataFrame(adata.layers['counts'].T) # get raw counts
    nUMI_spatial = pd.DataFrame(spatial_counts.sum(axis=0))

    # Convert coordinates and counts to R-friendly format
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_coords = ro.conversion.py2rpy(coords)  # Convert coordinates DataFrame to R
        r_spatial_counts = ro.conversion.py2rpy(pd.DataFrame(spatial_counts))  # Convert counts to R
        r_nUMI_spatial = ro.conversion.py2rpy(pd.Series(nUMI_spatial))  # Convert nUMI_spatial to R

    r_script = f"""
    library(spacexr)
    library(Matrix)
    
    reference <- loadRDS({ref_file})
    
    # Load spatial data
    
    puck <- SpatialRNA({r_coords}, {r_spatial_counts}, {r_nUMI_spatial})
    
    myRCTD <- create.RCTD(puck, reference, max_cores=4)
    myRCTD <- run.RCTD(myRCTD, doublet_mode='full_mode')
    
    saveRDS(myRCTD, "{sample_output_dir}/RCTD_results.rds")
    write.csv(myRCTD@results$weights, "{sample_output_dir}/RCTD_weights.csv")
    """
    
    ro.r(r_script)

# Load references from annotated single-cell data
single_cell_ref_h5 = "/storage/gge/Quique/TabulaeParalytica/single/GSE234774.h5"
single_cell_ref_h5_dir = "/storage/gge/Quique/TabulaeParalytica/single/"

adata_ref = sc.read(single_cell_ref_h5)
# Split adata into 2 based on injured (7d, 14d) and healthy individuals:

adata_ref_healthy = adata_ref[adata_ref.obs['label'].isin(['uninjured'])].copy()
adata_ref_injured = adata_ref[adata_ref.obs['label'].isin(['7d', '14d'])].copy()

# Make refs in appropriate format for RCTD
make_ref(adata_ref_healthy, single_cell_ref_h5_dir)
make_ref(adata_ref_injured, single_cell_ref_h5_dir)

for file_name in os.listdir(os.path.join(args.base_dir, "indiv_samples")):
    sample_path = os.path.join(args.base_dir, "indiv_samples", file_name, "outs", "matrices")
    for file in os.listdir(sample_path):
        if file.endswith('manual_delineation.h5ad'):
            sample_name = file_name
            print(f"Processing sample {sample_name}")
            # Define the directories to save plots and matrices
            sample_output_dir = os.path.join(output_dir, sample_name)
            os.makedirs(sample_output_dir, exist_ok=True)
            adata = sc.read_h5ad(os.path.join(sample_path, file))
            condition = adata.obs['condition'].unique()[0]
            deconvolution(adata, sample_name, condition, sample_output_dir)

