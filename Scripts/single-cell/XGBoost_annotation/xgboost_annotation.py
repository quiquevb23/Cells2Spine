'''
    Script to annotate cells based on the trained XGBoost model
'''

import scanpy as sc
import matplotlib.pyplot as plt
import os
import numpy as np
import argparse
import joblib
import seaborn as sns
import pandas as pd
from sklearn.feature_selection import mutual_info_classif
from sklearn.preprocessing import LabelEncoder
from xgboost import XGBClassifier
from sklearn.metrics import classification_report, confusion_matrix

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
model_dir = os.path.join(parent_dir, "XGBoost_models")
os.makedirs(model_dir, exist_ok=True)

output_base_dir = args.output_base_dir
parent_output_dir = os.path.dirname(output_base_dir)
output_dir = os.path.join(parent_output_dir, "XGBoost_outputs")
os.makedirs(output_dir, exist_ok=True)

single_cell_ref_h5 = "/storage/gge/Quique/TabulaeParalytica/single/GSE234774.h5"


def annotate_with_XGBoost(adata, single_cell_ref_h5, sample_name, model_dir, output_dir, sample_path):
    adata_ref = sc.read(single_cell_ref_h5)
    if sample_name == "Sample_1":    
        model_path = os.path.join(model_dir, "best_xgb_cell_type_model_healthy.joblib")
        model = joblib.load(model_path)
        adata_ref_subs = adata_ref[adata_ref.obs['label'].isin(['uninjured'])].copy()
    elif sample_name == "Sample_3":
        model_path = os.path.join(model_dir, "best_xgb_cell_type_model_injured.joblib")
        model = joblib.load(model_path)
        adata_ref_subs = adata_ref[adata_ref.obs['label'].isin(['7d', '14d'])].copy()
    
    adata.X = adata.raw.X # assign raw data
    # Normalize and log-transform again
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Step 4: Ensure the feature names are correctly set (this should have been done during training)
    # You can check that the model's booster has the correct feature names
    print("Feature names in model:", model.get_booster().feature_names)

    # Step 5: Select the same top-k genes from the new data (as done during training)
    # Get the top-k genes based on the feature importance
    top_k_genes = model.get_booster().feature_names  # These are the feature names from the trained model

    # Ensure the new data has the same genes in the same order as the training data
    X_new_selected = adata[:, top_k_genes].X  # Select the same top-k genes

    # Step 6: Ensure the new data is in a dense format (XGBoost requires dense format)
    X_new_dense = X_new_selected.toarray() if hasattr(X_new_selected, 'toarray') else X_new_selected

    # Step 7: Load the label encoder (used during training) to decode predictions
    label_encoder = LabelEncoder()
    label_encoder.fit(adata_ref_subs.obs['cell_type'])  # Fit it on the training data's labels

    # Predict cell annotations
    predictions = model.predict(X_new_dense)

    # Step 9: Map predicted numeric labels back to original labels
    y_pred_original = label_encoder.inverse_transform(predictions)

    # Step 10: Store or print the predictions
    adata.obs['predicted_annotations'] = y_pred_original

    # Print or save the predictions
    print(f"Predictions for new data: {y_pred_original}")

    # Plot UMAP of Leiden clusters and XGBoost predictions
    sc.tl.umap(adata)
    sc.pl.umap(adata, color = ['leiden_0.4', 'predicted_annotations'], wspace=0.5) #can use predicted_labels instead of majority voting
    plt.savefig(os.path.join(output_dir, f"XGBoost_UMAP_{sample_name}.png"), bbox_inches='tight')
    plt.close()

    df = pd.DataFrame({
        'reference': adata.obs['leiden_0.4'],
        'prediction': adata.obs['predicted_annotations']
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

    plt.savefig(os.path.join(output_dir, f"XGBoost_dotplot_{sample_name}.png"), bbox_inches='tight')
    plt.close()

    # Calculate genes differentially expressed for each group of predicted cell types
    sc.tl.rank_genes_groups(adata, groupby = 'predicted_annotations', method='wilcoxon') #this will show the DEA of cell-type groups assigned by celltypist
    sc.pl.rank_genes_groups_dotplot(adata, groupby = 'predicted_annotations', n_genes=5, standard_scale='var')
    plt.savefig(os.path.join(output_dir, f"XGBoost_rank_genes_groups_predicted_annotations_{sample_name}.png"), bbox_inches='tight')
    plt.close()

    # Create a crosstab to calculate the overlap
    crosstab = pd.crosstab(df['reference'], df['prediction'])
    crosstab_percent = crosstab.div(crosstab.sum(axis=1), axis=0) * 100
    
    # Plot the crosstab as a heatmap
    plt.figure(figsize=(10, 8))  # Adjust figure size as needed
    sns.heatmap(crosstab_percent, annot=True, fmt="d", cmap="viridis", cbar=True)

    # Add labels and title
    plt.xlabel('Predicted XGBoost Annotations')
    plt.ylabel('Leiden 0.4 Clusters')
    plt.title('Overlap Between Leiden 0.4 and Predicted Labels')

    plt.savefig(os.path.join(output_dir, f"XGBoost_overlap_heatmap_{sample_name}.png"), bbox_inches='tight')
    plt.close()


    # Save the updated AnnData object
    adata.write_h5ad(os.path.join(sample_path, f"{sample_name}_annotated_XGBoost.h5ad"))

    print(f"Annotations added and saved")


# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name, "outs", "matrices")
    for file in os.listdir(sample_path):
        if file.endswith('clustering.h5ad'):
            sample_name = file_name
            output_dir = os.path.join(output_base_dir, sample_name)

            adata = sc.read_h5ad(os.path.join(sample_path, file))
            annotate_with_XGBoost(adata, single_cell_ref_h5, sample_name, model_dir, output_dir, sample_path)
            #annotate(adata, sample_name, single_cell_ref_h5, model_dir, output_dir)