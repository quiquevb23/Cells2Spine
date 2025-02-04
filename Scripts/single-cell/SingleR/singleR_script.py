'''
    Script to perform annotation of single-cell data with SingleR using T.Paralytica reference
'''

import scanpy as sc
import argparse
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from scipy.sparse import issparse
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
parent_dir = os.path.dirname(base_dir)
model_dir = os.path.join(parent_dir, "celltypist_models")
os.makedirs(model_dir, exist_ok=True)

output_base_dir = args.output_base_dir

# Ensure that renv is activated in the R environment
ro.r('''
    library(renv)
    renv::restore()
''')

# Load R packages and functions
ro.r('''
    library(SingleR)
    library(SummarizedExperiment)
    run_singleR <- function(query_matrix, ref_matrix, ref_labels) {
        # Convert query and reference matrices to SummarizedExperiment
        query_se <- SummarizedExperiment(assays = list(counts = query_matrix))
        ref_se <- SummarizedExperiment(assays = list(counts = ref_matrix))
        
        # Run SingleR
        results <- SingleR(test = query_se, ref = ref_se, labels = ref_labels)
        
        # Return annotations
        return(data.frame(Cell = rownames(results), Labels = results$labels))
    }
''')

def annotate_with_singleR(adata, sample_name, single_cell_ref_h5, model_dir, output_dir, sample_path):
    """
    query_matrix: Pandas DataFrame (query single-cell dataset)
    ref_matrix: Pandas DataFrame (reference dataset)
    ref_labels: Pandas Series (cell type labels for reference)
    """
    adata_ref = sc.read(single_cell_ref_h5)
    if sample_name == "Sample_1":
        adata_ref_subs = adata_ref[adata_ref.obs['label'].isin(['uninjured'])].copy()
    elif sample_name == "Sample_3":
        adata_ref_subs = adata_ref[adata_ref.obs['label'].isin(['7d', '14d'])].copy()
    
    # Apply log1p normalization to adata_ref_subs
    # Adata query is already log1p normalized
    sc.pp.normalize_total(adata_ref_subs, target_sum=10000)
    sc.pp.log1p(adata_ref_subs)
    '''
        Adata ref contains annotated cell types at different resolution levels:
        - "cell_type": 15 cell types; "cell_l3": 22; "cell_l4": 44, "cell_l5": 94. 
    
    '''
    adata_ref_subs.var_names_make_unique()
    adata_ref_subs.obs_names_make_unique()
    # Step 1: Find the intersection of genes
    shared_genes = adata.var_names.intersection(adata_ref_subs.var_names)
    
    # Step 2: Subset the query and reference data
    query_subset = adata[:, shared_genes].copy()
    reference_subset = adata_ref_subs[:, shared_genes].copy()

    # Step 3: Convert data to Pandas DataFrame for R
    query_matrix = pd.DataFrame(query_subset.X.toarray(), index=query_subset.obs_names, columns=query_subset.var_names)
    ref_matrix = pd.DataFrame(reference_subset.X, index=reference_subset.obs_names, columns=reference_subset.var_names)
    
    query_matrix.index = query_matrix.index.str.upper()
    ref_matrix.index = ref_matrix.index.str.upper()

    common_genes = set(query_matrix.index).intersection(set(ref_matrix.index))
    
    # Subset both matrices
    query_matrix = query_matrix.loc[common_genes, :]
    ref_matrix = ref_matrix.loc[common_genes, :]

    assert query_matrix.shape[0] > 0, "Query matrix has no common genes."
    assert ref_matrix.shape[0] > 0, "Reference matrix has no common genes."

    # Extract cell type labels from the reference AnnData
    ref_labels = reference_subset.obs['cell_type'] 

    print("printing length common genes and reference labels")
    print(len(common_genes))
    print(len(ref_labels))
    # Step 4: Define and run the SingleR function in R
    r_code = """
    library(SingleR)
    library(SummarizedExperiment)

    run_singleR <- function(query_matrix, ref_matrix, ref_labels) {
        # Convert matrices to SummarizedExperiment objects
        query_se <- SummarizedExperiment(assays = list(logcounts = as.matrix(query_matrix)))
        ref_se <- SummarizedExperiment(assays = list(logcounts = as.matrix(ref_matrix)))
        
        # Run SingleR
        results <- SingleR(test = query_se, ref = ref_se, labels = ref_labels)
        
        # Return the annotations
        return(data.frame(Cell = rownames(results), Labels = results$labels))
    }
    """
    ro.r(r_code)
    run_singleR = ro.globalenv['run_singleR']

    # Convert data to R and run SingleR
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_query_matrix = ro.conversion.py2rpy(query_matrix)
        r_ref_matrix = ro.conversion.py2rpy(ref_matrix)
        r_ref_labels = ro.conversion.py2rpy(ref_labels)
    
    result_r = run_singleR(r_query_matrix, r_ref_matrix, r_ref_labels)
    
    with localconverter(ro.default_converter + pandas2ri.converter):
        # Convert the result back to a Pandas DataFrame
        result_df = ro.conversion.rpy2py(result_r)
    
    # Step 5: Add annotations to the full query AnnData object
    adata.obs['SingleR_Labels'] = result_df.set_index('Cell').reindex(adata.obs_names)['Labels']
    adata.write_h5ad(os.path.join(sample_path, f"{sample_name}_annotated_all.h5ad"))
    
    print(f"Completed annotation to {sample_name}")

single_cell_ref_h5 = "/storage/gge/Quique/TabulaeParalytica/single/GSE234774.h5"

# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name, "outs", "matrices")
    for file in os.listdir(sample_path):
        if file.endswith('both.h5ad'):
            sample_name = file_name
            output_dir = os.path.join(output_base_dir, sample_name, "SingleR")
            os.makedirs(output_dir, exist_ok=True)

            adata = sc.read_h5ad(os.path.join(sample_path, file))
            adata.var_names_make_unique()
            adata.obs_names_make_unique()
            annotate_with_singleR(adata, sample_name, single_cell_ref_h5, model_dir, output_dir, sample_path)
            #annotate(adata, sample_name, single_cell_ref_h5, model_dir, output_dir)