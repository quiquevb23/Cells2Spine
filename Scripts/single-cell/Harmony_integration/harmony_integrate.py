'''
	Script to integrate single cell count matrices from 2 batches

'''
import scanpy as sc
import anndata as ad
import pandas as pd
import os
import matplotlib.pyplot as plt
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
output_base_dir = args.output_base_dir

#Create new folders for "joined" datasets by Harmony
parent_dir = os.path.dirname(base_dir)
parent_output_dir = os.path.dirname(output_base_dir)

joined_base_dir = os.path.join(parent_dir, "joined")
joined_output_base_dir = os.path.join(parent_output_dir, "joined")

os.makedirs(joined_base_dir, exist_ok=True)
os.makedirs(joined_output_base_dir, exist_ok=True)

datas = []

for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name, "outs", "matrices")
    for file in os.listdir(sample_path):
        if file.endswith('PCAs.h5ad'):
            adata = sc.read_h5ad(os.path.join(sample_path, file))
            print(adata.shape)
            adata.X = adata.layers['counts'] # Restore raw counts for concatenation
            adata.obs['sample'] = file_name
            datas.append(adata)

# I need to return back to raw counts for concatenation, to correctly apply normalization
adata_combined = sc.concat(datas, index_unique = "_")
print(f"Size of adata_combined before filtering out cells with 0 counts: {adata_combined.shape}")
sc.pp.filter_cells(adata_combined, min_counts=1)
sc.pp.filter_genes(adata_combined, min_counts=1)
print(f"Size of adata_combined after filtering out cells with 0 counts: {adata_combined.shape}")

# After concatenation of samples, we need to preprocess data again, normalize and everything
# We still have raw counts stored in .layers['counts']
sc.pp.normalize_total(adata_combined) # The right way to normalize count data according to benchmark$
sc.pp.log1p(adata_combined)
# Select top 5000 HVG for combined dataset and filter them
sc.pp.highly_variable_genes(adata_combined, n_top_genes=5000, flavor="seurat")
sc.pl.highly_variable_genes(adata_combined)
plt.savefig(os.path.join(joined_output_base_dir, 'HVG_genes_combined.png'), bbox_inches='tight')
if adata_combined.raw is not None:
    print("Adata.raw exists")
else:
    print("Adata.raw not exists")
# We will not filter out HVG, but rather give them to input for PCA
sc.pp.scale(adata_combined) # Scale before PCA: all genes despite only using HVGs for PCA
sc.pp.pca(adata_combined, svd_solver="arpack", mask_var="highly_variable")
sc.pl.pca_scatter(adata_combined, color="sample")
plt.savefig(os.path.join(joined_output_base_dir, "PCA_scatter.png"))
plt.close()

sc.pl.pca_variance_ratio(adata_combined)
plt.savefig(os.path.join(joined_output_base_dir, "PCA_variance_ratio.png"))
plt.close()

print(adata_combined.shape)
print(adata_combined)

# We will compute now UMAP on combined dataset without integration
adata_noint = adata_combined.copy()
sc.pp.neighbors(adata_noint)
sc.tl.umap(adata_noint)
sc.pl.umap(adata_noint, color=["sample"], wspace=0.5)
plt.savefig(os.path.join(joined_output_base_dir, "UMAP_unintegrated.png"), bbox_inches='tight')
plt.close()
adata_noint.write(os.path.join(joined_base_dir, "adata_no_integrated.h5ad"))

# Now run Harmony
# Since we have normalized and scaled data after concat, we can run now Harmony, with first 50 PCs (as calculated)
sc.external.pp.harmony_integrate(adata_combined, key='sample')
print(adata_combined.obsm['X_pca_harmony'].shape)
sc.pp.neighbors(adata_combined, use_rep="X_pca_harmony")
sc.tl.umap(adata_combined)
sc.pl.umap(adata_combined, color=["sample"], wspace=0.5)
plt.savefig(os.path.join(joined_output_base_dir, "UMAP_integrated.png"), bbox_inches='tight')
plt.close()
adata_combined.write(os.path.join(joined_base_dir, "adata_integrated.h5ad"))
