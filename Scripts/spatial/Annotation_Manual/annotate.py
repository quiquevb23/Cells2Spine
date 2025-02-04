'''
	Script to annotate manually the spots for Visium spatial data, preprocessed with stSME
	It includes some layers
'''
import os
import pandas as pd
from sklearn import metrics
import multiprocessing as mp
import matplotlib.pyplot as plt
import argparse
import scanpy as sc
import scipy.sparse
import numpy as np
import seaborn as sns
from matplotlib import gridspec

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

joined_base_dir = os.path.join(parent_dir, "Manual_ann")
joined_output_base_dir = os.path.join(parent_output_dir, "Manual_ann")

os.makedirs(joined_base_dir, exist_ok=True)
os.makedirs(joined_output_base_dir, exist_ok=True)


marker_genes_1 = {
    "Neuron": ["Meg3"],
    "Astrocyte": ["Aqp4"],
    "OPC": ["Pdgfra"],
    "Microglia": ["C1qc"],
    "Oligodendrocyte": ["Cldn11"],
    "Endothelial": ["Vtn"],
}


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

marker_genes_special = {
    "Neuron": ["Atf3", "Gap43", "Sprr1a", "Neurod1"],  # Injury response in neurons
    "Astrocyte": ["Gfap", "Serpina3n", "Lcn2", "C3"],  # Reactive astrocytes after injury
    "OPC (Oligodendrocyte Precursor Cell)": ["Pdgfra", "Cspg4", "Sox10", "Vcan"],  # Reactive OPCs, scar-forming role
    "ODC (Mature Oligodendrocyte)": ["Mbp", "Mog", "Cnp", "Tcf7l2"],  # Oligodendrocyte demyelination/remyelination
    "Microglia": ["Cx3cr1", "Aif1", "Cd68", "Trem2"],  # Activated microglia/macrophage during injury
    "Endothelial": ["Flt1", "Cldn5", "Icam1", "Vegfa"],  # Blood-brain barrier disruption, angiogenesis
    "Pericyte": ["Pdgfrb", "Acta2", "Anpep", "Cspg4"],  # Scar formation and blood-brain barrier repair
    "Ependymal": ["Foxj1", "Vim", "Sox2", "Nestin"],  # Ependymal cell activation and neurogenesis potential
    "Stromal": ["Dcn", "Col1a1", "Tnc", "Pdgfrb"],  # Reactive stromal cells, involved in scar tissue
    "Leukocyte (Immune Cells)": ["Cd3e", "Ptprc", "Itgax", "Ccl2"],  # Infiltrating T cells, immune response
    "Neutrophil": ["S100a8", "S100a9", "Trem1", "Mpo"],  # Acute inflammation after injury
    "Macrophage": ["Cd68", "Cd163", "Mrc1", "Nos2"],  # Macrophage polarization and activation after SCI
}

'''
Neuron: Genes like Atf3, Gap43, and Sprr1a are typically upregulated in neurons after injury as part of the regenerative response and axon growth attempts.
Astrocyte: Reactive astrocytes after SCI show elevated levels of Gfap, Serpina3n, and Lcn2, which are involved in scar formation and immune modulation.
OPC: Pdgfra and Vcan are markers indicating the proliferation of oligodendrocyte precursor cells, often in an attempt to remyelinate after demyelination caused by injury.
ODC: Markers like Mbp and Tcf7l2 are related to remyelination and mature oligodendrocytes trying to repair lost myelin.
Microglia: Trem2 and Cd68 indicate microglial activation, playing a major role in inflammation and phagocytosis at the injury site.
Endothelial: Markers like Vegfa and Icam1 are associated with blood-brain barrier breakdown and subsequent angiogenesis post-injury.
Pericyte: These cells contribute to scar formation and express markers like Pdgfrb and Cspg4, responding to the vascular damage.
Ependymal: SCI induces activation of ependymal cells that express Nestin and Foxj1, with a potential role in regenerating the spinal cord tissue.
Fibroblast/Stromal: These cells, marked by Col1a1 and Tnc, participate in the formation of fibrotic scar tissue that blocks regeneration.
Leukocyte: T-cell infiltration (e.g., Cd3e) and other immune cells contribute to inflammation and modulate injury progression.
Neutrophil: S100a8 and Mpo are elevated early after injury, contributing to acute inflammation.
Macrophage: SCI drives the activation of macrophages, where markers like Cd68, Mrc1, and Nos2 indicate their pro-inflammatory or anti-inflammatory phenotypes.
'''

cell_types_of_interest = ["Ependymal", "Stromal", "Leukocyte (Immune Cells)", "Macrophage", "Microglia"]

for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name, "outs", "matrices")
    for file in os.listdir(sample_path):
        if file.startswith('adata_stlearn_domains_'): # we take the stSME norm counts
            sample_name = file_name
            adata = sc.read_h5ad(os.path.join(sample_path, file))

            # Set up the grid layout for 5 cell types in a single row
            fig = plt.figure(figsize=(25, 5))  # Adjust width (25) as needed
            gs = gridspec.GridSpec(1, 5, wspace=0.4)  # 1 row, 5 columns

            for idx, cell_type in enumerate(cell_types_of_interest):
                genes = marker_genes_special.get(cell_type, [])
                
                # Filter genes that exist in the adata
                existing_genes = [gene for gene in genes if gene in adata.var_names]

                if existing_genes:  # Only proceed if there are genes to plot
                    # Create a new column in adata for the aggregate expression of the cell type
                    adata.obs[cell_type] = adata.X[:, adata.var_names.isin(existing_genes)].sum(axis=1)

                    # Plot using the aggregate expression of the cell type
                    ax = fig.add_subplot(gs[idx])
                    sc.pl.spatial(
                        adata, 
                        color=cell_type,  # Use the new column we just created
                        title=f"{cell_type}", 
                        frameon=False, 
                        ax=ax,  # Plot in the correct subplot
                        show=False  # Avoid showing now, we'll save it later
                    )
                    # Display the genes as text in the plot
                    gene_text = ", ".join(existing_genes)
                    ax.text(0.5, -0.1, gene_text, fontsize='small', ha='center', transform=ax.transAxes, wrap=True)

            # Save the combined plot
            filename = f"combined_{sample_name}.png"
            filepath = os.path.join(joined_output_base_dir, filename)
            os.makedirs(os.path.dirname(filepath), exist_ok=True)
            plt.savefig(filepath)  # Save the figure
            plt.close()  # Close the plot to avoid display and free memory

            '''
            # Plotting all genes collectively
            for cell_type, genes in marker_genes_special.items():
                # Filter genes that exist in the adata
                existing_genes = [gene for gene in genes if gene in adata.var_names]
                
                if existing_genes:  # Only proceed if there are genes to plot
                    # Create a new column in adata for the aggregate expression of the cell type
                    adata.obs[cell_type] = adata.X[:, adata.var_names.isin(existing_genes)].sum(axis=1)

                    # Plot using the aggregate expression of the cell type
                    sc.pl.spatial(
                        adata, 
                        color=cell_type,  # Use the new column we just created
                        title=f"Spatial Expression of {cell_type} ({sample_name})", 
                        frameon=False, 
                        show=True
                    )
                    
                    # Save the figure
                    filename = f"{cell_type}_{sample_name}.png"
                    filepath = os.path.join(joined_output_base_dir, filename)
                    os.makedirs(os.path.dirname(filepath), exist_ok=True)
                    plt.savefig(filepath)  # Save the figure
                    plt.close()  # Close the plot to avoid display and free memory
            '''
