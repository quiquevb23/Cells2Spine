"""
    Script to compare annotations by Celltypist and Ingest and integrate them
"""


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

# Calculate cluster-cell type percentages for each method
def calculate_percentage_matrix(data):
    """
    Calculates the percentage matrix for cluster-cell type combinations.
    Args:
    - data: A DataFrame with columns 'cluster' and 'cell_type'.
    
    Returns:
    - percentage_matrix: A DataFrame with clusters as rows and cell types as columns, containing the percentage of cells assigned to each cell type in each cluster.
    """
    # Group by cluster and cell type and count occurrences
    cluster_cell_counts = data.groupby(['leiden_0.4', 'predicted_labels' if 'predicted_labels' in data.columns else 'cell_type']).size().unstack(fill_value=0)
    cluster_totals = cluster_cell_counts.sum(axis=1)  # Sum of cells in each cluster
    percentage_matrix = cluster_cell_counts.div(cluster_totals, axis=0) * 100  # Convert to percentages
    return percentage_matrix

def plot_agreement_heatmap(agreement_matrix, output_path):
    plt.figure(figsize=(10, 8))
    sns.heatmap(agreement_matrix, annot=True, cmap="YlGnBu", fmt=".1f", cbar_kws={"label": "Agreement (%)"})
    plt.title("Agreement Heatmap (Celltypist vs Ingest)")
    plt.xlabel("Cell Types")
    plt.ylabel("Clusters")
    plt.savefig(output_path)
    plt.close()  # Close the plot to avoid overlap with subsequent plots

def plot_disagreement_heatmap(disagreement_matrix, output_path):
    plt.figure(figsize=(10, 8))
    sns.heatmap(disagreement_matrix, annot=True, cmap="Reds", fmt=".1f", cbar_kws={"label": "Disagreement (%)"})
    plt.title("Disagreement Heatmap (Celltypist vs Ingest)")
    plt.xlabel("Cell Types")
    plt.ylabel("Clusters")
    plt.savefig(output_path)
    plt.close()  # Close the plot to avoid overlap with subsequent plots

def plot_scatter_plot(matrix_celltypist, matrix_ingest, output_path):
    # Ensure both matrices have the same clusters and cell types
    # Get the union of clusters and cell types
    all_clusters = matrix_celltypist.index.union(matrix_ingest.index)
    all_cell_types = matrix_celltypist.columns.union(matrix_ingest.columns)

    # Reindex both matrices to have the same clusters and cell types, filling missing values with 0
    matrix_celltypist_filled = matrix_celltypist.reindex(index=all_clusters, columns=all_cell_types, fill_value=0)
    matrix_ingest_filled = matrix_ingest.reindex(index=all_clusters, columns=all_cell_types, fill_value=0)

    # Flatten the matrices for scatter plot
    x = matrix_celltypist_filled.values.flatten()
    y = matrix_ingest_filled.values.flatten()

    # Scatter plot
    plt.figure(figsize=(8, 8))
    plt.scatter(x, y, alpha=0.7, c="blue")
    plt.plot([0, 100], [0, 100], "k--", label="Perfect Agreement")  # Diagonal line (perfect agreement)

    # Annotate points with large disagreement
    threshold = 20  # Adjust this threshold as needed
    for i, (xi, yi) in enumerate(zip(x, y)):
        if abs(xi - yi) > threshold:  # Highlight large disagreements
            plt.annotate(f"{i}", (xi, yi), fontsize=8, alpha=0.7)

    plt.xlabel("Celltypist Percentage")
    plt.ylabel("Ingest Percentage")
    plt.title("Scatter Plot: Celltypist vs Ingest Percentages")
    plt.legend()
    plt.savefig(output_path)
    plt.close()  # Close the plot to avoid overlap with subsequent plots

from scipy.stats import entropy

def calculate_entropy(matrix):
    """
    Calculate the entropy for each cluster based on the distribution of cell types.
    
    Args:
    - matrix (pd.DataFrame): A matrix of cluster-cell type percentages, where rows are clusters and columns are cell types.
    
    Returns:
    - pd.Series: A series with entropy values for each cluster.
    """
    # Normalize by the sum of each row to get proportions
    normalized_matrix = matrix.div(matrix.sum(axis=1), axis=0)
    
    # Calculate entropy for each row (cluster) based on the normalized proportions
    return normalized_matrix.apply(lambda x: entropy(x, base=2), axis=1)

def plot_entropy_plot(entropy_df, output_path):
    entropy_df.plot(kind="bar", figsize=(10, 6), colormap="viridis")
    plt.title("Cluster Entropy (Mixedness) for Each Annotation Method")
    plt.ylabel("Entropy")
    plt.xlabel("Clusters")
    plt.savefig(output_path)
    plt.close()  # Close the plot to avoid overlap with subsequent plots    

def compare_annotations(adata_celltypist, adata_ingest, sample_name, output_dir):
    """
    Compares annotations from two methods (Celltypist and Ingest) and generates percentage matrices.
    Args:
    - adata_celltypist: Annotated AnnData object from Celltypist method.
    - adata_ingest: Annotated AnnData object from Ingest method.
    - sample_name: The name of the sample for output file naming.
    - output_dir: Directory where the output files will be saved.
    """
    # Extract the relevant columns
    celltypist_clusters = adata_celltypist.obs[['leiden_0.4', 'predicted_labels']]
    ingest_clusters = adata_ingest.obs[['leiden_0.4', 'cell_type']]

    # Generate percentage matrices
    matrix_celltypist = calculate_percentage_matrix(celltypist_clusters)
    matrix_ingest = calculate_percentage_matrix(ingest_clusters)

    # Generate the agreement and disagreement matrices
    agreement_matrix = np.minimum(matrix_celltypist, matrix_ingest)
    disagreement_matrix = abs(matrix_celltypist - matrix_ingest)

    # Create output paths for saving the plots
    agreement_heatmap_path = os.path.join(output_dir, f"{sample_name}_agreement_heatmap.png")
    disagreement_heatmap_path = os.path.join(output_dir, f"{sample_name}_disagreement_heatmap.png")
    scatter_plot_path = os.path.join(output_dir, f"{sample_name}_scatter_plot.png")
    entropy_plot_path = os.path.join(output_dir, f"{sample_name}_entropy_plot.png")

    # Save the plots
    plot_agreement_heatmap(agreement_matrix, agreement_heatmap_path)
    plot_disagreement_heatmap(disagreement_matrix, disagreement_heatmap_path)
    plot_scatter_plot(matrix_celltypist, matrix_ingest, scatter_plot_path)
    
    # Calculate entropy and save the entropy plot
    entropy_celltypist = calculate_entropy(matrix_celltypist)
    entropy_ingest = calculate_entropy(matrix_ingest)
    entropy_df = pd.DataFrame({
        "Celltypist": entropy_celltypist,
        "Ingest": entropy_ingest
    })
    plot_entropy_plot(entropy_df, entropy_plot_path)

    print(f"Plots saved for {sample_name}:")
    print(f"- Agreement Heatmap: {agreement_heatmap_path}")
    print(f"- Disagreement Heatmap: {disagreement_heatmap_path}")
    print(f"- Scatter Plot: {scatter_plot_path}")
    print(f"- Entropy Plot: {entropy_plot_path}")

def add_annotations(adata_celltypist, adata_ingest, sample_path, file_name):
    adata_celltypist.obs['ingest_celltypes'] = adata_ingest.obs['cell_type']
    adata_celltypist.write_h5ad(os.path.join(sample_path, f"{file_name}_annotated_both.h5ad"))
    print("adata_saved")
# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name, "outs", "matrices")
    adata_celltypist = None
    adata_ingest = None
    output_dir = os.path.join(output_base_dir, file_name, "Comparison_Annotations")
    os.makedirs(output_dir, exist_ok=True)

    # Loop through the files in the sample directory
    for file in os.listdir(sample_path):
        if file.endswith('annotated.h5ad'):
            adata_celltypist = sc.read_h5ad(os.path.join(sample_path, file))
        elif file.endswith('annotated_ingest.h5ad'):
            adata_ingest = sc.read_h5ad(os.path.join(sample_path, file))
    add_annotations(adata_celltypist, adata_ingest, sample_path, file_name)

    # Ensure both Celltypist and Ingest annotations were loaded before comparison
    if adata_celltypist is not None and adata_ingest is not None:
        compare_annotations(adata_celltypist, adata_ingest, file_name, output_dir)
    else:
        print(f"Missing annotation files for sample: {file_name}")
