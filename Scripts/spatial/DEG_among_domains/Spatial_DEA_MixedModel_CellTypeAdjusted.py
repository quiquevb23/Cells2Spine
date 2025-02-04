
import os
import scanpy as sc
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
import argparse
import scipy.sparse
import stlearn as st
import seaborn as sns

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
output_base_dir = args.output_base_dir

#Create new folders for "joined" datasets by Harmony
parent_dir = os.path.dirname(base_dir)
parent_output_dir = os.path.dirname(output_base_dir)

joined_base_dir = os.path.join(parent_dir, "DEGs")
joined_output_base_dir = os.path.join(parent_output_dir, "DEGs")

os.makedirs(joined_base_dir, exist_ok=True)
os.makedirs(joined_output_base_dir, exist_ok=True)


def reorder_indices(adata, sample):
    int_sample = int(sample[-1]) #get sample name
    adata.obs.index = [f"{idx}_{int_sample}" for idx in adata.obs.index]
    return adata

def stSME_normalization(data, sample):
    tile_path = os.path.join("./tmp/tiles", sample)
    os.makedirs(tile_path, exist_ok=True)

    st.pp.filter_genes(data,min_cells=1)
    st.pp.normalize_total(data)
    st.pp.log1p(data)
    st.pp.tiling(data, tile_path)

    # this step uses deep learning model to extract high-level features from tile images
    # may need few minutes to be completed
    st.pp.extract_feature(data)
    # run PCA for gene expression data
    st.em.run_pca(data,n_comps=50)
    data_SME = data.copy()
    # apply stSME to normalise log transformed data
    st.spatial.SME.SME_normalize(data_SME, use_data="raw")
    data_SME.X = data_SME.obsm['raw_SME_normalized']

    return data_SME

def normalize(adata):
    sc.pp.normalize_total(adata, target_sum=1e4)
    # Logarithmize the data
    sc.pp.log1p(adata)
    # Scale the data (zero mean, unit variance)
    #sc.pp.scale(adata)
    return adata

def process(adata_combined):
    condition = adata_combined.obs['condition'].unique()[0]
    if condition == 'injured':
        subset_adata = adata_combined[adata_combined.obs['leiden_0.4'] == '6'].copy() #leiden cluste$
        subset_adata = normalize(subset_adata)
        return subset_adata
    else:
        subset_adata = adata_combined[adata_combined.obs['highlight_group'] == 'highlight'].copy() #copy selection for dorsal site of healthy
        subset_adata = normalize(subset_adata)
        return subset_adata

def estimate_coefficients(adata, gene_list, factor_names, expected_expression):
    """
    Estimate cell-type specific coefficients for each gene in each pixel and sample using a Negative Binomial GLM.
    
    Parameters:
    adata (AnnData): AnnData object containing the spatial counts.
    gene_list (list): List of genes to analyze.
    factor_names (list): List of cell type proportion column names in `adata.obs`.
    expected_expression (pd.DataFrame): DataFrame containing expected gene expression rates for each gene and cell type.
    
    Returns:
    results (dict): Dictionary with cell-type specific coefficients for each gene.
    """
    results = {}

    # Loop over each gene and fit the model
    for gene in gene_list:
        # Extract spatial counts for the gene
        counts = adata.raw[:, gene].X.toarray().flatten()  # Convert to 1D array
        
        # Initialize results for the gene
        gene_results = {}
        
        # Loop over each cell type
        for cell_type in factor_names:
            # Expected expression rate for the current gene and cell type
            proportions = adata.obs[cell_type].values

            # Create the design matrix with an intercept
            design_matrix = pd.DataFrame({
                'intercept': np.ones(len(counts)),
                'condition': adata.obs['condition'],
                'proportions': proportions
            })

            # Create dummy variables for the condition
            design_matrix = pd.get_dummies(design_matrix, columns=['condition'], drop_first=True)
            
            for col in design_matrix.columns:
                if "condition_" in col:
                    design_matrix[f'{col}_x_proportions'] = design_matrix[col] * design_matrix['proportions']
            design_matrix = design_matrix.astype(float)
    
            response = counts.astype(float)

            # Fit the GLM
            model = sm.GLM(
                response,
                design_matrix,
                family=sm.families.NegativeBinomial()
            )

            # Fit the model and handle potential errors
            try:
                results_gene_type = model.fit()
                # Store the coefficients for the current cell type
                gene_results[cell_type] = results_gene_type.params
            except Exception as e:
                print(f"Error fitting model for gene {gene} and cell type {cell_type}: {e}")

        results[gene] = gene_results
    
    return results



# Define the file path and sample names for cell type proportions
cell2loc_file_path = "/storage/gge/Quique/Cells2SpineData/Pilot/spatial/matrices/Deconvolution/cell2location_map"

datas_healthy = []
datas_injured = []
adata_dict = {}
datas_dorsal = []

# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name, "outs", "matrices")
    for file in os.listdir(sample_path):
        if file.endswith('_highlight.h5ad'):
            sample_name = file_name
            print(f"Processing healthy adata {sample_name}")
            # Define the directories to save plots and matrices
            output_dir = os.path.join(output_base_dir, sample_name)
            os.makedirs(output_dir, exist_ok=True)
            # This has the annotation for clusters but is not normalized

            adata = sc.read_h5ad(os.path.join(sample_path, file))
            sample = adata.obs['sample'].unique()[0]
            adata = reorder_indices(adata, sample)
            print(sample)
            # Map sample prefix to the correct spSpatial file
            if sample == "Spatial_1":
                spatial_file = "spSpatial_1"
            elif sample == "Spatial_2":
                spatial_file = "spSpatial_2"
            # Load cell type proportions from the mapped spSpatial file
            #adata = load_spSpatial_proportions(adata, cell2loc_file_path, spatial_file)
            #adata = st.convert_scanpy(adata)

            #adata_dict[sample_name] = stSME_normalization(adata, sample_name)
            adata_dict[sample_name] = adata
            datas_healthy.append(adata_dict[sample_name])

        elif file.startswith('adata_stlearn_common_domains_'):
            sample_name = file_name
            # Define the directories to save plots and matrices
            output_dir = os.path.join(output_base_dir, sample_name)
            os.makedirs(output_dir, exist_ok=True)
            # This has the annotation for clusters but is not normalized

            adata = sc.read_h5ad(os.path.join(sample_path, file))
            condition = adata.obs['condition'].unique()[0]
            if condition == 'injured':
                sample = adata.obs['sample'].unique()[0]
                adata = reorder_indices(adata, sample)
                print(sample)
                #adata = st.convert_scanpy(adata)
                print(f"Processing injured adata {sample_name}")

                # Map sample prefix to the correct spSpatial file
                if sample == "Spatial_3":
                    spatial_file = "spSpatial_3"
                elif sample == "Spatial_4":
                    spatial_file = "spSpatial_4"
                # Load cell type proportions from the mapped spSpatial file
                #adata = load_spSpatial_proportions(adata, cell2loc_file_path, spatial_file)

                #adata_dict[sample_name] = stSME_normalization(adata, sample_name)
                adata_dict[sample_name] = adata
                datas_injured.append(adata_dict[sample_name])


adata_healthy = sc.concat(datas_healthy, join='outer')
adata_injured = sc.concat(datas_injured, join='outer')

#Process to get only dorsal regions
healthy_subset = process(adata_healthy)
injured_subset = process(adata_injured)

# Assuming adata_healthy and adata_injured are AnnData objects with dorsal ROI subsetted.
# Merge data from both conditions into a single AnnData object
all_data = sc.concat([healthy_subset, injured_subset], join='outer')

#SC signatures
adata_ref = sc.read_h5ad('/storage/gge/Quique/Cells2SpineData/Pilot/spatial/matrices/Deconvolution/reference_signatures/sc.h5ad')

# Extract expected expression values from `adata_ref`
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' 
                                                         for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}' 
                              for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']

inf_aver.to_csv("./inf_aver.csv")

# Example gene list, cell type proportions, and expected expression
gene_list = list(all_data.var_names.intersection(inf_aver.index))

# Set up expected_expression DataFrame for coefficient estimation
expected_expression = inf_aver.loc[gene_list]

# Here add the celltype proportions
celltype_prop = pd.read_csv("./celltype_proportions.csv", index_col=0)

# Check if indices align, if not, handle mismatches as needed
if not celltype_prop.index.isin(all_data.obs.index).all():
    print("Warning: Some indices in `celltype_prop` are not found in `all_data.obs`.")

# Step 2: Define or extract factor names (you may need to set this manually if not available)
factor_names = celltype_prop.columns.tolist()  # Use columns from celltype_prop as factor names

if 'mod' not in all_data.uns:
    all_data.uns['mod'] = {}

# Step 5: Assign the factor names to all_data.uns['mod']['factor_names']
all_data.uns['mod']['factor_names'] = factor_names

# Add the cell type proportions into `all_data.obs` under the factor names
all_data.obs[all_data.uns['mod']['factor_names']] = celltype_prop[factor_names].reindex(all_data.obs.index)

all_data.write_h5ad("./alldata_proportions.h5ad")

# Run the coefficient estimation
results = estimate_coefficients(all_data, gene_list, factor_names, expected_expression)

# Flatten the results dictionary into a DataFrame
flattened_results = []

for gene, cell_types in results.items():
    for cell_type, coefficients in cell_types.items():
        row = {
            'gene': gene,
            'cell_type': cell_type,
            **coefficients  # Unpack the coefficients for each gene and cell type
        }
        flattened_results.append(row)

# Convert the list of dictionaries to a DataFrame
results_df = pd.DataFrame(flattened_results)

# Save the DataFrame to a CSV file
results_df.to_csv(os.path.join(joined_output_base_dir,"cell_type_specific_DEG_coefficients.csv"), index=False)


'''
# Step 1: Prepare DataFrame with gene expression, cell type proportions, and condition
# Extract gene expression data as a DataFrame
expression_df = pd.DataFrame(all_data.X.toarray(), columns=all_data.var_names, index=all_data.obs.index)

# Add cell type proportions and metadata (e.g., condition and replicate)
expression_df['condition'] = all_data.obs['condition']
expression_df['sample'] = all_data.obs['sample_id']  # assuming 'sample_id' indicates replicates
cell_type_proportions = pd.DataFrame(all_data.obsm['cell_type_proportions'], index=all_data.obs.index)
expression_df = pd.concat([expression_df, cell_type_proportions], axis=1)

# Step 2: Define and Fit the Mixed Linear Model for Each Gene
# Iterate through each gene and fit a mixed model
results = {}
for gene in all_data.var_names:
    formula = f"{gene} ~ condition + " + " + ".join(cell_type_proportions.columns) + " + (1|sample)"
    model = smf.mixedlm(formula, expression_df, groups=expression_df["sample"])
    result = model.fit(reml=False)  # Set REML to False for maximum likelihood estimation
    results[gene] = {
        'p_value': result.pvalues['condition[T.injured]'],
        'logFC': result.params['condition[T.injured]'],
        'std_error': result.bse['condition[T.injured]']
    }

# Convert results to a DataFrame for easy analysis and filtering
results_df = pd.DataFrame(results).T  # Transpose to have genes as rows
results_df['adj_p_value'] = sm.stats.multipletests(results_df['p_value'], method='fdr_bh')[1]

# Step 3: Filter DEGs Based on Adjusted p-value and logFC thresholds
degs = results_df[(results_df['adj_p_value'] < 0.05) & (np.abs(results_df['logFC']) > 1)]

# Step 4: Map DEGs back to Cell Type Contributions
degs['cell_type_contribution'] = degs['logFC'] * cell_type_proportions.mean(axis=0)

print(degs[['logFC', 'adj_p_value', 'cell_type_contribution']])
'''
