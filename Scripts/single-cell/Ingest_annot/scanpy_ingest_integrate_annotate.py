import scanpy as sc
import matplotlib.pyplot as plt
import os
import pandas as pd
import numpy as np
import argparse
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

'''
BASE_DIR='/storage/gge/Quique/Cells2SpineData/Pilot/single-cell/full_seq_aggr/indiv_samples'
OUTPUT_BASE_DIR='/home/quiquevb/Cells2Spine/Cells2Spine/Outputs/single-cell/full_seq_aggr/indiv_samples'
'''

base_dir = args.base_dir
parent_dir = os.path.dirname(base_dir)

output_base_dir = args.output_base_dir

# Path for reference of single cell data
single_cell_ref_h5 = "/storage/gge/Quique/TabulaeParalytica/single/GSE234774.h5"

def normalize(adata):
    # Adata_ref_subs is raw
    sc.pp.normalize_total(adata, target_sum=10000)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=5000)
    sc.pp.scale(adata)
    sc.pp.pca(adata, svd_solver="arpack", mask_var="highly_variable")
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    return adata

def run_ingest(adata, single_cell_ref_h5, sample_name, output_dir, sample_path):
    adata_ref = sc.read(single_cell_ref_h5)
    if sample_name == "Sample_1":    
        adata_ref_subs = adata_ref[adata_ref.obs['label'].isin(['uninjured'])].copy()
    elif sample_name == "Sample_3":
        adata_ref_subs = adata_ref[adata_ref.obs['label'].isin(['7d', '14d'])].copy()

    if adata_ref_subs.raw is not None:
        print("adata_ref.raw exists")
    else:
        print("adata_ref.raw does not exist")
    adata_ref_subs.var_names_make_unique()
    common_vars = adata.var_names.intersection(adata_ref_subs.var_names)
    adata.X = adata.layers['counts'].copy() # Back to raw
    adata = adata[:, common_vars].copy()
    adata_ref_subs = adata_ref_subs[:, common_vars].copy()
    # We have adata processed with neighborhood embeddings and UMAP 
    adata = normalize(adata)
    adata_ref_subs = normalize(adata_ref_subs)
    
    sc.tl.ingest(adata, adata_ref_subs, obs='cell_type')

    # Plot UMAP of Leiden clusters and XGBoost predictions
    sc.pl.umap(adata, color = ['leiden_0.4', 'cell_type'], wspace=0.5) #can use predicted_labels instead of majority voting
    plt.savefig(os.path.join(output_dir, f"Ingest_UMAP_{sample_name}.png"), bbox_inches='tight')
    plt.close()

    df = pd.DataFrame({
        'reference': adata.obs['leiden_0.4'],
        'prediction': adata.obs['cell_type']
    })
    count_table = df.groupby(['reference', 'prediction']).size().reset_index(name='cell_count')
    
    # Normalize to get the fraction of cells in each reference group
    reference_totals = count_table.groupby('reference')['cell_count'].transform('sum')
    count_table['fraction'] = count_table['cell_count'] / reference_totals

    # Pivot the table for plotting
    dotplot_data = count_table.pivot(index='prediction', columns='reference', values='fraction').fillna(0)

    # Plot the dotplot
    plt.figure(figsize=(10, 6))
    sns.heatmap(dotplot_data, annot=True, fmt=".2f", cmap="Blues", cbar_kws={'label': 'Fraction of cells'})
    plt.title("Fraction of Cells by Reference and Predictions")
    plt.xlabel("Reference Groups (Leiden Clusters)")
    plt.ylabel("Predicted Cell Types")
    plt.tight_layout()

    plt.savefig(os.path.join(output_dir, f"Ingest_dotplot_{sample_name}.png"), bbox_inches='tight')
    plt.close()

    group_sizes = adata.obs['cell_type'].value_counts()
    # Reassign rare groups to "Rare" category
    adata.obs['cell_type'] = adata.obs['cell_type'].apply(
        lambda x: x if group_sizes[x] > 1 else 'Rare'
    )

    # Differential expression analysis: Filter out groups with fewer than 2 cells
    def rank_genes_with_check(adata, groupby, output_path):
        # Get counts per group
        group_counts = adata.obs[groupby].value_counts()

        # Filter groups with fewer than 2 cells
        valid_groups = group_counts[group_counts > 1].index.tolist()

        # Subset adata to only valid groups
        adata_valid = adata[adata.obs[groupby].isin(valid_groups)]

       	if len(valid_groups) > 1:
            # Perform DEA
            sc.tl.rank_genes_groups(adata_valid, groupby=groupby, method='wilcoxon')
            sc.pl.rank_genes_groups_dotplot(adata_valid, groupby=groupby, n_genes=5, standard_scale='var')
            plt.savefig(output_path, bbox_inches='tight')
            plt.close()
        else:
            print(f"Skipping DEA for {groupby} due to insufficient group sizes (only one group with multiple cells)")

    # Perform DEA for predicted labels
    rank_genes_with_check(adata, 'cell_type', os.path.join(output_dir, f"Ingest_rank_genes_groups_predicted_annotations_{sample_name}.png"))

    # Create a crosstab to calculate the overlap
    crosstab = pd.crosstab(df['reference'], df['prediction'])
    crosstab_percent = crosstab.div(crosstab.sum(axis=1), axis=0) * 100
    
    # Plot the crosstab as a heatmap
    plt.figure(figsize=(10, 8))  # Adjust figure size as needed
    sns.heatmap(crosstab_percent, annot=True, cmap="viridis", cbar=True)

    # Add labels and title
    plt.xlabel('Predicted Ingest Annotations')
    plt.ylabel('Leiden 0.4 Clusters')
    plt.title('Overlap Between Leiden 0.4 and Predicted Labels')

    plt.savefig(os.path.join(output_dir, f"Ingest_overlap_heatmap_{sample_name}.png"), bbox_inches='tight')
    plt.close()


    # Save the updated AnnData object
    adata.write_h5ad(os.path.join(sample_path, f"{sample_name}_annotated_ingest.h5ad"))



# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name, "outs", "matrices")
    for file in os.listdir(sample_path):
        if file.endswith('clustering.h5ad'):
            sample_name = file_name
            
            output_dir = os.path.join(output_base_dir, sample_name, "Ingest")
            os.makedirs(output_dir, exist_ok=True)
            adata = sc.read_h5ad(os.path.join(sample_path, file))
            
            run_ingest(adata, single_cell_ref_h5, sample_name, output_dir, sample_path)
            #annotate(adata, sample_name, single_cell_ref_h5, model_dir, output_dir)
