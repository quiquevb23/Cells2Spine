'''
	Script to integrate single cell count matrices from 2 batches

'''
import scanpy as sc
import anndata as ad
import pandas as pd
import os
import matplotlib.pyplot as plt

base_dir = '/storage/gge/Quique/Cells2SpineData/Pilot/single-cell/matrices/'
output_base_dir = '/home/quiquevb/Cells2Spine/Cells2Spine/Outputs/single-cell/joined'

adata_list = []

for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name)
    for file in os.listdir(sample_path):
        if file.endswith('normalized.h5ad'):
            file_path = os.path.join(sample_path, file)
            adata = sc.read_h5ad(file_path)
            adata.obs['sample'] = file_name
            adata_list.append(adata)

adata_combined = adata_list[0].concat(*adata_list[1:], batch_key="sample", join="outer")

sc.pp.pca(adata_combined)
adata_noint = adata_combined.copy()
sc.pp.neighbors(adata_noint, use_rep='X_pca')
sc.tl.umap(adata_noint)
sc.pl.umap(adata_noint, color=[sample], wspace=1)
plt.savefig(os.path.join(output_base_dir, "UMAP_unintegrated.png"))
adata_noint.write(os.path.join(output_base_dir, "adata_no_integrated.h5ad"))

#Now run Harmony
sc.external.pp.harmony_integrate(adata_combined, key='sample')
sc.pp.neighbors(adata_combined, use_rep="X_pca_harmony")
sc.tl.umap(adata_combined)
sc.pl.umap(adata_combined, color=[sample], wspace=1)
plt.savefig(os.path.join(output_base_dir, "UMAP_integrated.png"))
adata_combined.write(os.path.join(output_base_dir, "adata_integrated.h5ad"))

