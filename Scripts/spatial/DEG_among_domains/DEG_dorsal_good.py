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

base_dir = os.path.join(args.base_dir, "DEA_dorsal")
output_dir = os.path.join(args.output_base_dir, "DEA_dorsal")

os.makedirs(base_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

# Define paths for files for adata ST data deconvoluted either with separated sc references or single sc

adata_files = '/storage/gge/Quique/Cells2SpineData/Pilot/spatial/matrices/indiv_samples/Spatial_1/outs/matrices'


def add_deconv_scores(adata, sample_name, sample_output_dir):
    path = "/storage/gge/Quique/Cells2SpineData/Pilot/spatial/matrices/Deconvolution/cell2location_map"
    for file in os.listdir(path):
        if file.endswith(f"{sample_name}_ct_clusters.h5ad"):
            deconv_file = sc.read_h5ad(os.path.join(path, file)) # get adata file with deconvolution scores
    # 
    #check if indices match first
    indices1 = adata.obs.index
    indices2 = deconv_file.obs.index

    # Check if indices match
    if indices1.equals(indices2):
        print("The indices match!")
        
        # Add deconvolution scores
        cell_type_columns = deconv_file.uns['mod']['factor_names']

        #Assign cell-type values
        adata.obs[cell_type_columns] = deconv_file.obs[cell_type_columns]
        print("Successfully assigned the selected columns from deconv_file.obs to adata.obs.")
        
        # Add cell-type clusters
        adata.obs['ct_leiden_0.4'] = deconv_file.obs['ct_leiden_0.4']
        adata.obs['ct_leiden_0.8'] = deconv_file.obs['ct_leiden_0.8']
        
        if sample_name in ['Spatial_3', 'Spatial_4']:
            # Extract macrophage density values (replace 'macrophage_density_column' with the actual column name)
            if 'Macrophages, chemotaxis-inducing' in adata.obs.columns:
                # Extract macrophage density values
                macrophage_density = adata.obs['Macrophages, chemotaxis-inducing']
            else:
                print('Macrophage, chemotaxis-inducing do not exist')
                return

            print(adata.obs['Macrophages, chemotaxis-inducing'].describe())
            print(adata.obs['Macrophages, chemotaxis-inducing'].head())

            # Calculate the 90th and 95th quantile thresholds
            q90 = macrophage_density.quantile(0.90)
            q95 = macrophage_density.quantile(0.95)
            print(f"90th Quantile: {q90}")
            print(f"95th Quantile: {q95}")

            # Identify macrophage-rich groups for each threshold
            adata.obs['macrophage_rich_q90'] = macrophage_density > q90
            adata.obs['macrophage_rich_q95'] = macrophage_density > q95
            
            # Plot spatial data for the 90th quantile group
            sc.pl.spatial(
                adata, 
                color='macrophage_rich_q90', 
                title=f"Macrophage-rich (90th Percentile): {sample_name}"
            )
            plt.savefig(os.path.join(sample_output_dir, f'macrophage_rich_q90_{sample_name}.png'))
            plt.close()

            # Plot spatial data for the 95th quantile group
            sc.pl.spatial(
                adata, 
                color='macrophage_rich_q95', 
                title=f"Macrophage-rich (95th Percentile): {sample_name}"
            )           
            plt.savefig(os.path.join(sample_output_dir, f'macrophage_rich_q95_{sample_name}.png'))
            plt.close()
 
        return adata
    else:
        print("The indices do not match.")

def select_highlight_group(adata, sample_name, sample_output_dir):
    spatial_coords = adata.obsm['spatial']
    # Find min and max for x (column 0) and y (column 1)
    x_min = np.min(spatial_coords[:, 0])
    x_max = np.max(spatial_coords[:, 0])
    y_min = np.min(spatial_coords[:, 1])
    y_max = np.max(spatial_coords[:, 1])

    # Spot diameter provided: 118.07816
    if sample_name == 'Spatial_1':
        spot_diameter = 118.07816
    elif sample_name == 'Spatial_2':
        spot_diameter = 118.0915099
    elif sample_name == 'Spatial_3':
        spot_diameter = 118.12852000000001
    elif sample_name == 'Spatial_4':
        spot_diameter = 118.11344000000001
   

    # We want grids of 2x2 spots
    grouping_factor = 1

    # Calculate step size for each grid based on the diameter of the spots
    x_grid_step = spot_diameter * grouping_factor
    y_grid_step = spot_diameter * grouping_factor

    # Create grid IDs by dividing x and y coordinates by the step size and flooring the result
    grid_x_id = np.floor((spatial_coords[:, 0] - np.min(spatial_coords[:, 0])) / x_grid_step).astype(int)
    grid_y_id = np.floor((spatial_coords[:, 1] - np.min(spatial_coords[:, 1])) / y_grid_step).astype(int)
    
    # Combine the x and y grid IDs into a single group ID
    adata.obs['grid_group'] = (grid_x_id * 1000 + grid_y_id).astype(str)  # Assign unique IDs
    
    if sample_name == 'Spatial_1':
        #grids for adata_1
        grids_1 = ['16003', '16005', '17003', '17004', '17006', '18003', '18005', '18007', '18008', '18010', '18011',
                    '20003', '20004', '20006', '20007', '20009', '20010', '20012', '20013', '21003', '21005', '21007',
                    '21008', '21010', '21011', '21013', '21014', '22009', '22010', '22012', '22013'
                    ]
        # Highlight only the '3003' group, and set others to gray
        adata.obs['highlight_group'] = np.where(adata.obs['grid_group'].isin(grids_1), 'highlight', 'other')

    elif sample_name == 'Spatial_2':
        # grids for adata_2
        grids_2 = ['4023', '4025', '4027', '5023', '5024', '5026', '5027', '6023', '6025', '6027', '8024', '8026',
                '9023', '9025', '9027', '10023', '10024', '10026', '12023', '12025', '13023', '13024', '14022', 
                '14023', '16023', '16024', '17023'
                ]
        # Highlight only the '3003' group, and set others to gray
        adata.obs['highlight_group'] = np.where(adata.obs['grid_group'].isin(grids_2), 'highlight', 'other')
    
    elif sample_name == 'Spatial_3':
        # grids for adata_2
        grids_3 = ['21016', '21017',
                '22015', '22016', '22018',
                '24013', '24014', '24016', 
                '25012', '25013', '25015', '25016',
                '26011', '26013', '26014', '26016',
                '28010', '28012', '28015',
                '29010', '29011', '29013', '29014',
                '30007', '30009', '30010', '30012',
                '32006', '32008', '32010', '32011', 
                '33006', '33007', '33009',
                '34010']
        # Highlight only the '3003' group, and set others to gray
        adata.obs['highlight_group'] = np.where(adata.obs['grid_group'].isin(grids_3), 'highlight', 'other')
    
    elif sample_name == 'Spatial_4':
        # grids for adata_4
        grids_4 = ['16006', '16010', '16011',
                '17003', '17004', '17006', '17007', '17009', '17010', '17012', '17013',
                '18002', '18003', '18005', '18007', '18008', '18010', '18011', '18013', '18014',
                '20003', '20004', '20006', '20007', '20009', '20010','20012',
                '21002', '21003', '21005', '21007', '21008', '21010',
                '22003', '22004'
                ]
        # Highlight only the '3003' group, and set others to gray
        adata.obs['highlight_group'] = np.where(adata.obs['grid_group'].isin(grids_4), 'highlight', 'other')
    
    # Visualize the spots, highlighting the '3003' group
    sc.pl.spatial(adata, color=['highlight_group'], img_key='hires', 
              title='Highlighting Grid Group 3003',   # 'other' in gray, 'highlight' in red
              show=False)
    plt.savefig(os.path.join(sample_output_dir, 'Dorsal_highlight.png'))
    plt.close()

    return adata

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


# Create list adatas
adatas = []
# Import stlearn datasets with common clusters and all genes
for file_name in os.listdir(os.path.join(args.base_dir, "indiv_samples")):
    sample_path = os.path.join(args.base_dir, "indiv_samples", file_name, "outs", "matrices")
    for file in os.listdir(sample_path):
        if file.startswith('adata_stlearn_common'):
            sample_name = file_name
            print(f"Processing sample {sample_name}")
            # Define the directories to save plots and matrices
            sample_output_dir = os.path.join(output_dir, sample_name)
            os.makedirs(sample_output_dir, exist_ok=True)
            
            adata = sc.read_h5ad(os.path.join(sample_path, file))
            condition = adata.obs['condition'].unique()[0]
            # Add deconvolution scores based on matched indices
            adata = add_deconv_scores(adata, sample_name, sample_output_dir)

            adata = select_highlight_group(adata, sample_name, sample_output_dir)
            
            adata.write_h5ad(os.path.join(args.base_dir, "indiv_samples", file_name, "outs", "matrices", f"adata_{file_name}_dorsal_highlight.h5ad"))
            # Slice adata to select highlight
            if 'highlight_group' in adata.obs:
                adata_subset = adata[adata.obs['highlight_group'] == 'highlight'].copy()
                if adata_subset.raw is not None:
                    adata_subset.X = adata_subset.raw.X.copy() # Use raw counts for DEA
                else:
                    raise ValueError(f"Raw counts are not available in the AnnData object for sample {sample_name}.")
                # Append subset adata to list
                adatas.append(adata_subset)
            else:
                raise KeyError("The 'highlight_group' attribute is not present in `adata.obs`.")
            
            # Check spots that are outliers according to PCA in dorsal
            #check_outliers(adata_subset, sample_name, sample_output_dir)
            
# Combine all subsetted AnnData objects for pseudobulk aggregation
combined_adata = ad.concat(adatas)
combined_adata.obs_names_make_unique()

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

'''
# Filter out genes with zero variance across samples
non_variable_genes = grouped.var(axis=1) == 0
print(f"Number of non-variable genes: {non_variable_genes.sum()}")

grouped_filtered = grouped.loc[~non_variable_genes]

print("Grouped shape after removing non-variable genes:", grouped_filtered.shape)

# Check if any pseudobulk sample has all zeros
zero_counts_samples = (grouped_filtered == 0).sum(axis=0)  # Sum of zeros in each aggregated sample (column)
print("Zero counts in pseudobulk samples:")
print(zero_counts_samples[zero_counts_samples > 0])  # Prints samples with zero counts
'''

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

'''
# Assign row names (gene names) to DGEList counts matrix
rownames_r_counts = list(grouped_filtered.index)  # Ensure this is the correct gene name list

gene_names = list(rownames_r_counts)
print(len(gene_names)) #up to this point it works
'''
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

'''
# Filter the gene names using the same 'keep' vector
rfiltered_gene_names = [gene_names[i] for i in range(len(keep)) if keep[i]]

print(len(filtered_gene_names)) #this length should be ~4000genes


# Create a new DGEList object with the filtered counts and gene names
dge_filtered = edgeR.DGEList(counts=dge_filtered_counts, group=r_condition)

# Extract the counts matrix from the filtered DGEList object
dge_counts = dge_filtered.rx2("counts")

# Convert the counts matrix (from R) into a pandas DataFrame
with localconverter(ro.default_converter + pandas2ri.converter):
    dge_counts_df = ro.conversion.rpy2py(dge_counts)

# Check the shape of the DataFrame
print(f"Shape of dge_counts_df: {dge_counts_df.shape}")

# Convert to pandas DataFrame (if not already)
dge_counts_df = pd.DataFrame(dge_counts_df)

dge_counts_df.index = filtered_gene_names  # Add gene names from filtered gene list (excluding lowly ex$
dge_counts_df.reset_index(inplace=True)
dge_counts_df.rename(columns={"index": "gene"}, inplace=True)

# Save raw counts to CSV
raw_counts_csv_path = os.path.join(output_dir, "raw_counts.csv")
dge_counts_df.to_csv(raw_counts_csv_path, index=True)  # The index should correspond to gene names
print(f"Raw counts saved to {raw_counts_csv_path}")


# Estimate normalization factors
dge_filtered = edgeR.calcNormFactors(dge_filtered)

# Access the raw counts and normalization factors
raw_counts = dge_filtered.rx2("counts")
# Extract normalization factors
norm_factors = dge_filtered.rx2("samples").rx2("norm.factors")

# Print normalization factors
print(norm_factors)

normalized_counts = edgeR.cpm(dge_filtered, normalized_lib_sizes = True)
with localconverter(ro.default_converter + pandas2ri.converter):
    normalized_counts_df = ro.conversion.rpy2py(normalized_counts)

# Check the shape of the DataFrame
print(f"Shape of dge_counts_df: {normalized_counts_df.shape}")

# Convert to pandas DataFrame (if not already)
normalized_counts_df = pd.DataFrame(normalized_counts_df)

normalized_counts_df.index = filtered_gene_names  # Add gene names from filtered gene list (excluding lowly ex$
normalized_counts_df.reset_index(inplace=True)
normalized_counts_df.rename(columns={"index": "gene"}, inplace=True)

# Save raw counts to CSV
normalized_counts_df_csv_path = os.path.join(output_dir, "normalized_counts.csv")
normalized_counts_df.to_csv(normalized_counts_df_csv_path, index=True)  # The index should correspond to gene names
print(f"Raw counts saved to {normalized_counts_df_csv_path}")
'''
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
csv_output_path = os.path.join(output_dir, "DEA_dorsal_pseudobulk.csv")
results_df.to_csv(csv_output_path, index=False)
print(f"Results saved to {csv_output_path}")
