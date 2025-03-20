"""
	Subset dorsal region of adata for DEA among healthy and injured using pseudobulk approach
    1st: add information about cell2location from csv files
    2nd: performed matched dorsal DEA
"""

import os
import pandas as pd
from sklearn import metrics
import multiprocessing as mp
import matplotlib.pyplot as plt
import argparse
import scanpy as sc
import scipy.sparse
import numpy as np
import seaborn as sns
from rpy2 import robjects as ro
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
import anndata as ad

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

#base_dir = os.path.join(args.base_dir, "DEA_dorsal")
output_dir = os.path.join(args.output_base_dir, "DEA_areas")

#os.makedirs(base_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

# Define paths for files for adata ST data deconvoluted either with separated sc references or single sc

adata_files = '/storage/gge/Quique/Cells2SpineData/Pilot/spatial/matrices/indiv_samples/Spatial_1/outs/matrices'


def check_outliers(adata, sample_name, sample_output_dir):
    def process_and_pca(data, label):
        """Helper function to preprocess and perform PCA"""
        
        sc.pp.scale(data)  # Z-score scaling

        # Step 2: Perform PCA
        sc.tl.pca(data, svd_solver='arpack')
        pca_results = data.obsm['X_pca']

        # Step 3: Visualize PCA
        plt.figure(figsize=(8, 6))
        plt.scatter(pca_results[:, 0], pca_results[:, 1], c='blue', alpha=0.6, label=f'Spots ({label})')
        plt.title(f'PCA of Highlight Group ({label} Data)')
        plt.xlabel('PC1')
        plt.ylabel('PC2')
        plt.legend()
        plt.grid(True)

        # Save the plot
        output_path = f"{sample_output_dir}/{sample_name}_PCA_outliers_scaling_{label}.png"
        plt.savefig(output_path)
        plt.show()
        print(f"PCA plot saved for {label} data at {output_path}")

    # Check with Raw Counts
    if adata.raw is not None:
        print("Processing raw counts...")
        raw_data = adata.raw.to_adata()  # Create an AnnData object from raw counts
        process_and_pca(raw_data, "Raw_Counts")
    else:
        print("Warning: No raw counts found in adata.raw.")

    # Check with SME-Normalized Counts
    if "raw_SME_normalized" in adata.obsm:
        print("Processing SME-normalized counts...")
        sme_data = adata.copy()  # Use a copy of the AnnData object with SME-normalized counts
        process_and_pca(sme_data, "SME-Normalized")
    else:
        print("Warning: No SME-normalized data found in adata.obsm.")

def do_DEA(combined_adata, area_output_dir):
    # Aggregate counts to pseudobulk
    df = combined_adata.to_df()  # Extract counts matrix
    metadata = combined_adata.obs[["sample", "condition"]]  # Extract metadata

    if df.shape[0] < df.shape[1]:  # Check if it's (samples, genes) instead of (genes, samples)
        print("Transposing count matrix to match (genes, samples) format.")
        df = df.T  # Transpose to have genes as rows and samples as columns
    
    # Check initial shape and verify the data
    print("Initial df shape:", df.shape)  # Should be (genes, samples)
    print("First few rows of df:")
    print(df.head())
    
    # Check the metadata shape
    print("Metadata shape:", metadata.shape)  # Should be (samples, 2)

    # Print first few rows of metadata to ensure it's correct
    print(metadata.head())

    # Aggregate spots in each region and sample
    grouped = df.T.groupby(metadata['sample']).sum().T  # Transpose before groupby, then transpose back

    grouped_transposed = grouped.T
    # Transpose the grouped DataFrame so that the samples are rows (instead of columns)
    #grouped_transposed = grouped_filtered.T  # Transpose to get samples as rows

    # Create pseudobulk AnnData object
    pseudobulk_metadata = pd.DataFrame(grouped_transposed.index.tolist(), columns=["sample"])
    pseudobulk_metadata["condition"] = metadata.groupby("sample")["condition"].first().values

    # Check the result of the metadata
    print(pseudobulk_metadata.head())

    # Assign 'condition' to pseudobulk_metadata for the 4 aggregated samples
    pseudobulk_adata = ad.AnnData(
        X=grouped_transposed.values, 
        obs=pseudobulk_metadata,
        #var=pd.DataFrame(index=grouped_filtered.index)
        var=pd.DataFrame(index=grouped.index)
    )

    # Extract counts and metadata from pseudobulk AnnData
    counts = pd.DataFrame(
        pseudobulk_adata.X, 
        index=pseudobulk_adata.obs.index, 
        columns=pseudobulk_adata.var_names
    )

    metadata = pseudobulk_adata.obs
    print("counts shape and columns:")
    # Check the shape of the counts DataFrame before conversion
    print(counts.shape)  # It should be (128, 4) if there are 128 genes and 4 pseudobulk samples

    # Check column names of the counts DataFrame
    print(counts.columns)  # It should match the 4 pseudobulk samples

    # Simply transpose counts
    counts = counts.T
    print("Pseudobulk obs head:")
    print(pseudobulk_adata.obs.head())  # Check that the "condition" column is properly added

    # Activate the pandas-to-R DataFrame converter
    pandas2ri.activate()

    # Import edgeR and base R packages
    base = importr("base")
    utils = importr("utils")
    edgeR = importr("edgeR")

    # Convert Python DataFrame to R DataFrame
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_counts = ro.conversion.py2rpy(counts)
        r_metadata = ro.conversion.py2rpy(pseudobulk_metadata)
        r_genes = ro.conversion.py2rpy(pseudobulk_adata.var_names)

    # Ensure 'condition' is a factor in R
    r_condition = r_metadata.rx2('condition')  # Extract condition as a factor

    print("Sample conditions:", r_condition)  # Verify which condition corresponds to each sample

    # Create design matrix for the model
    design_matrix = ro.r('model.matrix')(ro.Formula('~ condition'), data=r_metadata)
    print("design_matrix")
    print(design_matrix)

    condition_levels = r_metadata.rx2('condition').levels
    print("Condition levels:", condition_levels)

    # Create DGEList object in R
    #dge = edgeR.DGEList(counts=r_counts, group=r_condition)
    dge = edgeR.DGEList(counts=r_counts, group=r_condition, genes=r_genes)
    # The following might also work
    #dge = edgeR.DGEList(counts=r_counts, group=design_matrix, genes=r_genes)

    # Print the counts before filtering
    dge_counts_before = r['as.matrix'](dge.rx2('counts'))
    print("Counts before filtering:")
    print(dge_counts_before)

    # Filter lowly expressed genes
    #keep = edgeR.filterByExpr(dge, group=r_condition, keep_lib_sizes=False)
    keep = edgeR.filterByExpr(dge, group=r_condition)

    # Print the number of genes kept
    print(f"Number of genes retained after filtering: {sum(keep)}")
    print(f"Total number of genes before filtering: {len(keep)}")

    # Assuming `dge` and `keep` are already defined
    r.assign('dge', dge)
    r.assign('keep', keep)

    # Subset the DGEList using R syntax
    dge_filtered_counts = r('dge[keep, , keep.lib.sizes=FALSE]')

    # Print the counts after filtering
    dge_counts_after = r['as.matrix'](dge_filtered_counts.rx2('counts'))
    print("Counts after filtering:")
    print(dge_counts_after)

    dge_filtered = edgeR.normLibSizes(dge_filtered_counts) #this already normalizes

    # Estimate dispersions
    ##dge_filtered = edgeR.estimateDisp(dge_filtered, tagwise=True, design=design_matrix)
    dge_filtered = edgeR.estimateDisp(dge_filtered)

    # Perform exact test for two-group comparison
    de_results = edgeR.exactTest(dge_filtered)

    # Extract results table
    results_table = edgeR.topTags(de_results, n=ro.r("nrow")(dge_filtered)).rx2("table")

    # Convert R DataFrame to Python DataFrame
    with localconverter(ro.default_converter + pandas2ri.converter):
        results_df = ro.conversion.rpy2py(results_table)

    '''
    # Ensure gene names are included in the output CSV
    results_df.index = filtered_gene_names  # Add gene names from filtered gene list (excluding lowly expressed)
    results_df.reset_index(inplace=True)
    results_df.rename(columns={"index": "gene"}, inplace=True)
    '''

    # Inspect the results
    print(results_df.head())

    # Save the results to CSV
    csv_output_path = os.path.join(area_output_dir, "DEA_dorsal_pseudobulk.csv")
    results_df.to_csv(csv_output_path, index=False)
    print(f"Results saved to {csv_output_path}")

# Create list adatas
adatas = []
# Import stlearn datasets with common clusters and all genes
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
            adatas.append(adata)

all_areas = set()
for adata in adatas:
    all_areas.update(adata.obs['manual_delineation'].unique())

for area in all_areas:
    # Create a list to hold subsets for each area in all adatas    
    area_adatas = []
    
    for adata in adatas:
        # Create a subset of adata for each area and add it to the list
        area_subset = adata[adata.obs['manual_delineation'].astype(str) == str(area)].copy()
        
        # Make sure you're using the raw counts from the 'counts' layer
        area_subset.X = area_subset.layers['counts'].copy()
        
        # Append the subset to the area_adatas list
        area_adatas.append(area_subset)

    combined_adata = ad.concat(area_adatas)
    combined_adata.obs_names_make_unique()

    # Define output dir for each area
    area_output_dir = os.path.join(output_dir, area)
    os.makedirs(area_output_dir, exist_ok=True)

    # Perform DEA
    do_DEA(combined_adata, area_output_dir)


