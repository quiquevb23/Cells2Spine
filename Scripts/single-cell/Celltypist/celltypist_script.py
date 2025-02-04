import scanpy as sc
import celltypist
from celltypist import models
import argparse
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os

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

def annotate_with_tabulae(adata, sample_name, single_cell_ref_h5, model_dir, output_dir, sample_path):
    adata_ref = sc.read(single_cell_ref_h5)
    if sample_name == "Sample_1":
        adata_ref_subs = adata_ref[adata_ref.obs['label'].isin(['uninjured'])].copy()
    elif sample_name == "Sample_3":
        adata_ref_subs = adata_ref[adata_ref.obs['label'].isin(['7d', '14d'])].copy()
    
    # Apply log1p normalization to adata_ref_subs
    sc.pp.normalize_total(adata_ref_subs, target_sum=1e4)
    sc.pp.log1p(adata_ref_subs)
    '''
        Adata ref contains annotated cell types at different resolution levels:
        - "cell_type": 15 cell types; "cell_l3": 22; "cell_l4": 44, "cell_l5": 94. 
    
    '''
    # Train the model for cell_type signature identification (understand how it works)
    new_model = celltypist.train(
        adata_ref_subs, 
        labels = 'cell_type', #choose between different levels of resolution 
        n_jobs = 10, 
        feature_selection = True
    )
    # Path for saving the model
    model_path = os.path.join(model_dir, f"model_for_{sample_name}.pkl")
    # First save the model locally
    new_model.write(model_path)

    adata.X = adata.layers['counts'].copy()
    # Apply log1p normalization to adata_ref_subs
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Now that model is saved, we specify the model path (recommended as this ensures the model is intact every time it is loaded).
    predictions = celltypist.annotate(adata, model = model_path, majority_voting = True)    
    adata = predictions.to_adata()
    sc.tl.umap(adata)
    sc.pl.umap(adata, color = ['leiden_0.4', 'predicted_labels', 'majority_voting'], wspace=0.5) #can use predicted_labels instead of majority voting
    plt.savefig(os.path.join(output_dir, f"celltypist_umap_{sample_name}.png"), bbox_inches='tight')
    plt.close()

    celltypist.dotplot(predictions, use_as_reference = 'leiden_0.4', use_as_prediction = 'majority_voting')
    plt.savefig(os.path.join(output_dir, f"celltypist_dotplot_{sample_name}.png"), bbox_inches='tight')
    plt.close()

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
    rank_genes_with_check(adata, 'predicted_labels', os.path.join(output_dir, f"rank_genes_groups_predicted_labels_{sample_name}.png"))

    # Perform DEA for leiden clusters
    rank_genes_with_check(adata, 'leiden_0.4', os.path.join(output_dir, f"rank_genes_groups_leiden_{sample_name}.png"))

    df = pd.DataFrame({
        'reference': adata.obs['leiden_0.4'],
        'prediction': adata.obs['predicted_labels']
    })

    # Create a crosstab to calculate the overlap
    crosstab = pd.crosstab(df['reference'], df['prediction'])
    crosstab_percent = crosstab.div(crosstab.sum(axis=1), axis=0) * 100

    count_table = df.groupby(['reference', 'prediction']).size().reset_index(name='cell_count')
    
    # Normalize to get the fraction of cells in each reference group
    reference_totals = count_table.groupby('reference')['cell_count'].transform('sum')
    count_table['fraction'] = count_table['cell_count'] / reference_totals

    # Plot the crosstab as a heatmap
    plt.figure(figsize=(10, 8))  # Adjust figure size as needed
    sns.heatmap(crosstab_percent, annot=True, cmap="viridis", cbar=True)

    # Add labels and title
    plt.xlabel('Predicted Ingest Annotations')
    plt.ylabel('Leiden 0.4 Clusters')
    plt.title('Overlap Between Leiden 0.4 and Predicted Labels')

    plt.savefig(os.path.join(output_dir, f"overlap_heatmap_{sample_name}.png"), bbox_inches='tight')
    plt.close()
    
    adata.write_h5ad(os.path.join(sample_path, f"{sample_name}_annotated.h5ad"))
    '''
    # Now we can check the celltype-driving genes in both training and query datasets
    model = models.Model.load(model_path)
    cell_types = model.cell_types
    for cell in cell_types:
        top_3_genes = model.extract_top_markers(cell, 3)
        print(top_3_genes)
        # Check expression of the three genes in the training set.
        sc.pl.violin(adata_ref_subs, top_3_genes, groupby = 'cell_type', rotation = 90)
        plt.savefig(os.path.join(output_dir, f"top3genes_trainingset_{cell}_{sample_name}.png"))
        # Check expression of the three genes in the query set.
        # Here we use `majority_voting` from CellTypist as the cell type labels for this dataset.
        sc.pl.violin(adata, top_3_genes, groupby = 'majority_voting', rotation = 90)
        plt.savefig(os.path.join(output_dir, f"top3genes_queryset_{cell}_{sample_name}.png"))
    '''

#def annotate(adata, sample_name, single_cell_ref_h5, model_dir, output_dir):


single_cell_ref_h5 = "/storage/gge/Quique/TabulaeParalytica/single/GSE234774.h5"

# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name, "outs", "matrices")
    for file in os.listdir(sample_path):
        if file.endswith('clustering.h5ad'):
            sample_name = file_name
            output_dir = os.path.join(output_base_dir, sample_name, "Celltypist")
            os.makedirs(output_dir, exist_ok=True)

            adata = sc.read_h5ad(os.path.join(sample_path, file))
            annotate_with_tabulae(adata, sample_name, single_cell_ref_h5, model_dir, output_dir, sample_path)
            #annotate(adata, sample_name, single_cell_ref_h5, model_dir, output_dir)
