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

# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name, "outs", "matrices")
    for file in os.listdir(sample_path):
        if file.endswith('clustering.h5ad'):
            sample_name = file_name
            output_dir = os.path.join(output_base_dir, sample_name)

            adata = sc.read_h5ad(os.path.join(sample_path, file))
            updated_marker_genes = {}

            # Loop through each cell type in the marker_genes dictionary
            for cell_type, genes in marker_genes.items():
                # Filter the genes that are present in adata.var_names
                present_genes = [gene for gene in genes if gene in adata.var_names]
    
                # If any genes are present, update the dictionary
                if present_genes:
                    updated_marker_genes[cell_type] = present_genes

            sc.pl.dotplot(
                adata,
                groupby="leiden_0.2",
                var_names=updated_marker_genes,
                standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
                show=False
            )
            plt.savefig(os.path.join(output_dir, 'dotplot_marker_genes_0.2.png'))
            plt.close()

            sc.pl.dotplot(
                adata,
                groupby="leiden_0.4",
                var_names=updated_marker_genes,
                standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
                show=False
            )
            plt.savefig(os.path.join(output_dir, 'dotplot_marker_genes_0.4.png'))
            plt.close()

            sc.pl.dotplot(
                adata,
                groupby="leiden_0.8",
                var_names=updated_marker_genes,
                standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
                show=False
            )
            plt.savefig(os.path.join(output_dir, 'dotplot_marker_genes_0.8.png'))
            plt.close()

            sc.tl.rank_genes_groups(adata, groupby="leiden_0.2", method="wilcoxon", use_raw=False)
            sc.pl.rank_genes_groups_dotplot(
                adata, groupby="leiden_0.2", standard_scale="var", n_genes=5
            )
            plt.savefig(os.path.join(output_dir, 'dotplot_DEGs_0.2.png'))
            plt.close()

            sc.tl.rank_genes_groups(adata, groupby="leiden_0.4", method="wilcoxon", use_raw=False)
            sc.pl.rank_genes_groups_dotplot(
                adata, groupby="leiden_0.4", standard_scale="var", n_genes=5
            )
            plt.savefig(os.path.join(output_dir, 'dotplot_DEGs_0.4.png'))
            plt.close()

            sc.tl.rank_genes_groups(adata, groupby="leiden_0.8", method="wilcoxon", use_raw=False)
            sc.pl.rank_genes_groups_dotplot(
                adata, groupby="leiden_0.8", standard_scale="var", n_genes=5
            )
            plt.savefig(os.path.join(output_dir, 'dotplot_DEGs_0.8.png'))
            plt.close()

            #If we find as DEGs genes such as Mito, Hb or Ribo, we may need to regress them out
            adata.write(os.path.join(sample_path, sample_name + "DEGs.h5ad"))

