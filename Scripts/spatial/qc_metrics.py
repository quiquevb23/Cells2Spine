import os
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import squidpy

# Define directories
base_dir = '/storage/gge/Quique/Cells2SpineData/Pilot/spatial/matrices'
output_base_dir = '/home/quiquevb/Cells2Spine/Cells2Spine/Outputs/spatial'

# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name)
    outs_path = os.path.join(sample_path, "outs")
    for file in os.listdir(outs_path):
        if file.endswith('filtered_feature_bc_matrix.h5'):
            sample_name = file.replace('_filtered_feature_bc_matrix.h5', '')
            counts_file = os.path.join(outs_path, file)
            
            # Create output directories
            output_dir = os.path.join(output_base_dir, sample_name, 'QC_metrics')
            os.makedirs(output_dir, exist_ok=True)
        
            # Load the data
            adata = squidpy.read.visium(path=outs_path, counts_file=counts_file)
            adata.var_names_make_unique()
            adata.layers['counts'] = adata.X
            # Basic QC metrics
            adata.var['mt'] = adata.var_names.str.startswith('mt-')  # Identify mitochondrial genes
            adata.var['ribo'] = adata.var_names.str.startswith(("Rps", "Rpl"))
            adata.var["hb"] = adata.var_names.str.contains(("Hbb", "Hba", "Hbq"))
            adata.raw = adata.copy()
            sc.pp.calculate_qc_metrics(
                adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
            )
        
            # Plot QC metrics
            p1 = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
            p1.savefig(os.path.join(output_dir, 'total_counts_displot.png'))
            #sc.pl.violin(adata, 'total_counts')
            sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],
              jitter=0.4, multi_panel=True, show=False)
            plt.savefig(os.path.join(output_dir, 'violin.png'))
            sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False)
            plt.savefig(os.path.join(output_dir, 'scatter.png'))

            #plot spatially the QC metrics
            plt.rcParams["figure.figsize"] = (8, 8)
            sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"], show=False)
            plt.savefig(os.path.join(output_dir, 'spatial_counts.png'))
            sc.pl.spatial(adata, img_key="hires", color=["pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"], show=False)
            plt.savefig(os.path.join(output_dir, 'spatial_pct.png'))
            #save adata
            adata.write(os.path.join(sample_path, sample_name + "qc_metrics.h5ad"))
        else:
            print("filtered_feature_matrix not exists")
