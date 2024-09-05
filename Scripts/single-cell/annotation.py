'''
        Script to normalize individual samples of single-cell and perform PCA and clustering
'''

import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import os
import scipy as sp

base_dir = '/storage/gge/Quique/Cells2SpineData/Pilot/single-cell/matrices/'
output_base_dir = '/home/quiquevb/Cells2Spine/Cells2Spine/Outputs/single-cell'

marker_genes = {
    "Neuron": ["Snap25", "Map2", "Rbfox3", "Syp"],
    "Astrocyte": ["Ntsr2", "Htra1", "Aqp4"],
    "OPC": ["Plp1", "Mobp", "Mag", "Mog"],
    "ODC": ["Gpr17", "Pdgfra", "Sox10"],
    "Microglia": ["Ctss", "Cx3cr1", "Aif1", "Ly86"],
    "Endothelial": ["Cldn5", "Flt1", "Tek", "Cd34","Pecam1", "Prom1"],
    "Pericyte": ["Pdgfrb", "Vtn", "Myl9"],
    "Ependyma": ["Foxj1", "Sox2", "Rsph1", "Ak7"],
    "Stromal": ["Dcn", "Apod", "Gsn", "Col1a1", "Col3a1"],
    "Erythrocyte": ["Hbb-bt", "Hba-a1", "Hba-a2"],
    "Leukocyte": ["Ms4a4b", "Ltb", "Ctsw", "Cd3e"], 
    "Neutrophil": ["S100a8", "S100a9", "Trem1"],
}

def filter_marker_genes(adata, marker_genes):
    filtered_marker_genes = {category: [gene for gene in genes if gene in adata.var_names]
                             for category, genes in marker_genes.items()}
    
    # Remove any categories that end up with no genes after filtering
    filtered_marker_genes = {category: genes for category, genes in filtered_marker_genes.items() if genes}
    
    return filtered_marker_genes    



# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name)
    for file in os.listdir(sample_path):
        if file.endswith('_clustering.h5ad'):
            output_dir = os.path.join(output_base_dir, sample_name, 'c')
            sample_name = file.replace('_clustering.h5ad', '')
            adata = sc.read_h5ad(os.path.join(sample_path, file))
            sc.pl.dotplot(
                adata_pp,
                groupby="leiden_0.05",
                var_names=marker_genes,
                standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
                show=False
            )
            plt.savefig(os.path.join(output_dir, 'dotplot_marker_genes_0.05.png'))
            sc.tl.rank_genes_groups(adata, groupby="leiden_0.05", method="wilcoxon")
            sc.pl.rank_genes_groups_dotplot(
                adata, groupby="leiden_0.05", standard_scale="var", n_genes=5
            )
            plt.savefig(os.path.join(output_dir, 'dotplot_DEGs_0.05.png'))
            #If we find as DEGs genes such as Mito, Hb or Ribo, we may need to regress them out
