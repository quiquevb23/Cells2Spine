import os
import torch
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp
from GraphST import GraphST
import matplotlib.pyplot as plt
import seaborn as sns

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

single_cell_ref_h5 = "/storage/gge/Quique/TabulaeParalytica/single/GSE234774.h5"

os.environ['R_HOME'] = '/home/quiquevb/.conda/envs/graphst/lib/R'



def deconvolution(adata, sample, output_dir):
    # preprocessing for ST data
    GraphST.preprocess(adata)
    # use stSME normalized data

    # build graph
    GraphST.construct_interaction(adata)
    GraphST.add_contrastive_label(adata)

    # read scRNA daa
    adata_sc = sc.read(single_cell_ref_h5)
    adata_sc.var_names_make_unique()
    # preprocessing for scRNA data
    GraphST.preprocess(adata_sc)
    
    # find overlap genes
    from GraphST.preprocess import filter_with_overlap_gene
    adata, adata_sc = filter_with_overlap_gene(adata, adata_sc)

    # get features
    GraphST.get_feature(adata)

    import torch
    # Run device, by default, the package is implemented on 'cpu'. We recommend using GPU.
    device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')

    # Train model
    model = GraphST.GraphST(adata, adata_sc, epochs=1200, random_seed=50, device=device, deconvolution=True)
    adata, adata_sc = model.train_map()

    # Project cells into spatial space
    from GraphST.utils import project_cell_to_spot
    # the function seems to not work correctly
    print(adata.obsm['map_matrix'])
    project_cell_to_spot(adata, adata_sc) # we will remove retain percent just in case
    print(adata)
    
    # Visualization of spatial distribution of scRNA-seq data
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    with mpl.rc_context({'axes.facecolor':  'black',
                         'figure.figsize': [4.5, 5]}):

            sc.pl.spatial(adata, cmap='magma',
                      # selected cell types
                      color=['Ventral', 'Dorsal', 'Microglia', 'Astrocytes', 'Peripheral immune cells'],
                      ncols=5, size=1.3,
                      img_key='hires',
                      # limit color scale at 99.2% quantile of cell abundance
                      vmin=0, vmax='p99.2',
                      show=False
                     )

            plt.savefig(os.path.join(output_dir, "deconvolution_{sample}.png"), dpi=300, bbox_inches='tight')
            plt.close()

    return adata


adata_dict = {}
# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name, "outs", "matrices")
    for file in os.listdir(sample_path):
        if file.startswith('adata_stlearn_domains'):
            sample = file_name
            # Define the directories to save plots and matrices
            output_dir = os.path.join(output_base_dir, file_name)
            os.makedirs(output_dir, exist_ok=True)
            # populate the dictionary with adatas
            adata_dict[sample] = sc.read_h5ad(os.path.join(sample_path, file))
            print(f"Processing sample: {sample}")
            adata_dict[sample] = deconvolution(adata_dict[sample], sample, output_dir)

