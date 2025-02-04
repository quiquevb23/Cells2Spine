'''
        Scrip to perform DEG among clusters and find GO enrichment
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
import stlearn as st
import gseapy as gp

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

joined_base_dir = os.path.join(parent_dir, "DEGs")
joined_output_base_dir = os.path.join(parent_output_dir, "DEGs")

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


def stSME_normalization(data, sample):
    tile_path = os.path.join("./tmp/tiles", sample)
    os.makedirs(tile_path, exist_ok=True)

    st.pp.filter_genes(data,min_cells=1)
    st.pp.normalize_total(data)
    st.pp.log1p(data)
    st.pp.tiling(data, tile_path)

    # this step uses deep learning model to extract high-level features from tile images
    # may need few minutes to be completed
    st.pp.extract_feature(data)
    # run PCA for gene expression data
    st.em.run_pca(data,n_comps=50)
    data_SME = data.copy()
    # apply stSME to normalise log transformed data
    st.spatial.SME.SME_normalize(data_SME, use_data="raw")
    data_SME.X = data_SME.obsm['raw_SME_normalized']

    return data_SME

def normalize(adata):
    sc.pp.normalize_total(adata, target_sum=1e4)
    # Logarithmize the data
    sc.pp.log1p(adata)
    # Scale the data (zero mean, unit variance)
    sc.pp.scale(adata)
    return adata

def process(adata_combined):
    condition = adata_combined.obs['condition'].unique()[0]
    if condition == 'injured':
        subset_adata = adata_combined[adata_combined.obs['leiden_0.4'] == '6'].copy() #leiden clusters are str
#       subset_adata = normalize(subset_adata)
        return subset_adata
    else:
        subset_adata = adata_combined[adata_combined.obs['leiden_0.4'].isin(['0', '1'])]
        subset_adata = normalize(subset_adata)
        return subset_adata

def DEG_injured(adata_filtered):
    print("Processing injured subset")
#    sc.pp.normalize_total(adata_filtered, inplace=True)

#    sc.pp.scale(adata_filtered)  # Optional: scale your data if not already done

    # Calculate the mean expression of each gene
    mean_expression = np.ravel(adata_filtered.X.mean(axis=0))  # Convert sparse matrix to array

    # Create a DataFrame for genes and their mean expression
    gene_means = pd.DataFrame({'gene': adata_filtered.var_names, 'mean_expression': mean_expression})

    # Sort genes by mean expression (or use log fold change, etc.)
    ranked_genes = gene_means.sort_values(by='mean_expression', ascending=False)

    # Prepare ranked gene list (gene: mean_expression)
    ranked_gene_list = ranked_genes.set_index('gene')['mean_expression']
    glist=ranked_gene_list.index[:1000].to_list()
#    print(glist)
    # Perform GSEA with the ranked gene list
    '''
    gsea_go_results = gp.prerank(
        rnk=ranked_gene_list,  # Ranked gene list
        gene_sets='GO_Biological_Process_2021',  # Change as needed
        outdir=joined_output_base_dir,  # Output directory
        min_size=3,  # Minimum size of gene sets
        max_size=2000,  # Maximum size of gene sets
        permutation_num=2000,  # Number of permutations
        format='png',  # Output format for plots
        seed=42  # Set seed for reproducibility
    )
    '''
    enr = gp.enrichr(glist, # or "./tests/data/gene_list.txt",
                 gene_sets=['KEGG_2019_Mouse'],
                 organism='Mouse', # don't forget to set organism to the one you desired! e.g. Yeast
                 outdir=joined_output_base_dir, # don't write to disk
                 cutoff=0.5
                 )
    print(enr.results.head(5))
    # simple plotting function
    from gseapy import barplot, dotplot
    enr_kegg_results_df = enr.results
    enr_kegg_results_df.to_csv(os.path.join(joined_output_base_dir, 'kegg_enrichment_results.csv'))  # S$

    # categorical scatterplot
    ax = dotplot(enr.results,
              column="Adjusted P-value",
              x='Gene_set', # set x axis, so you could do a multi-sample/library comparsion
              size=5,
              top_term=10,
              figsize=(3,5),
              title = "KEGG",
              xticklabels_rot=45, # rotate xtick labels
              cutoff=0.5
              )
    ax.figure.savefig(os.path.join(joined_output_base_dir, "kegg_enrichment_dotplot.png"), bbox_inches='tight')

    '''
    # Save GO results to a DataFrame and CSV
    go_results_df = gsea_go_results.res2d
    go_results_df.to_csv(os.path.join(joined_output_base_dir,'enrichment_results.csv'))  # Save results to CSV
    print(gsea_go_results.res2d.head(20))

    # Plot and save enrichment results for the top 10 GO terms
    top_terms_go = go_results_df.head(10)['Term']
    gp.plot.dotplot(gsea_go_results.res2d, column='FDR q-val', title='Top 10 GO Enriched Terms', 
                    cmap=plt.cm.viridis, figsize=(4,5), size=6, cutoff=0.25)
    plt.savefig(os.path.join(joined_output_base_dir, 'top_10_enriched_go_terms.png'), bbox_inches='tight')
    plt.close()
    gene_sets_list = gp.get_library_name()
    #print(gene_sets_list)
    # Perform GSEA with the ranked gene list for KEGG
    gsea_kegg_results = gp.prerank(
        rnk=ranked_gene_list,  # Ranked gene list
        gene_sets='KEGG_2019_Mouse',  # Change as needed (ensure this matches the organism)
        outdir=joined_output_base_dir,  # Output directory for KEGG results
        min_size=3,  # Minimum size of gene sets
        max_size=2000,  # Maximum size of gene sets
        permutation_num=2000,  # Number of permutations
        format='png',  # Output format for plots
        seed=42  # Set seed for reproducibility
    )

    # Save KEGG results to a DataFrame and CSV
    kegg_results_df = gsea_kegg_results.res2d
    kegg_results_df.to_csv(os.path.join(joined_output_base_dir, 'kegg_enrichment_results.csv'))  # Save results to CSV
    print(gsea_kegg_results.res2d.head(20))

    # Plot and save enrichment results for the top 10 KEGG terms
    top_terms_kegg = kegg_results_df.head(10)['Term']
    gp.plot.dotplot(gsea_kegg_results.res2d, column='FDR q-val', title='Top 10 KEGG Enriched Terms', 
                    cmap=plt.cm.viridis, figsize=(4,5), size=6, cutoff=0.25)
    plt.savefig(os.path.join(joined_output_base_dir, 'top_10_enriched_kegg_terms.png'), bbox_inches='tight')
    plt.close()
    '''

datas_healthy = []
datas_injured = []
adata_dict = {}

# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name, "outs", "matrices")
    for file in os.listdir(sample_path):
        if file.startswith('adata_stlearn_common_domains_'):
            sample_name = file_name

            # Define the directories to save plots and matrices
            output_dir = os.path.join(output_base_dir, sample_name)
            os.makedirs(output_dir, exist_ok=True)
            # This has the annotation for clusters but is not normalized
            
            adata = sc.read_h5ad(os.path.join(sample_path, file))
            adata = st.convert_scanpy(adata)
            adata_dict[sample_name] = stSME_normalization(adata, sample_name)
#            adata_dict[sample_name] = adata
            condition = adata.obs['condition'].unique()[0]
            if condition == 'healthy':
                datas_healthy.append(adata_dict[sample_name])
            elif condition == 'injured':
                datas_injured.append(adata_dict[sample_name])
            else:
                print("No condition found")
#adata_healthy = sc.concat(datas_healthy, join='inner', index_unique = "_")
adata_injured = sc.concat(datas_injured, join='inner', index_unique = "_")

#healthy_subset = process(adata_healthy)
injured_subset = process(adata_injured)

DEG_injured(injured_subset)

'''
# Combine both subsets
adata_combined_subset = sc.concat([injured_subset, healthy_subset])

sc.pl.violin(adata_combined_subset, ['n_genes_by_counts', 'total_counts'], groupby='condition', jitter=0.4)
plt.savefig(os.path.join(joined_output_base_dir, f"Counts_Genes_distribution.png"), bbox_inches='tight')
plt.close()

# Perform differential expression analysis between the two groups
sc.tl.rank_genes_groups(adata_combined_subset, groupby='condition', groups=['injured'], reference='healthy', method='t-test')

# View the results
sc.pl.rank_genes_groups(adata_combined_subset, n_genes=20, sharey=False)
plt.savefig(os.path.join(joined_output_base_dir, f"Rank_genes_groups.png"), bbox_inches='tight')
plt.close()

# Extract the results for the injured group
deg_results = adata_combined_subset.uns['rank_genes_groups']

# Top genes for injured group
top_genes_injured = deg_results['names']['injured']
print(top_genes_injured[:20])  # Top 100 DE genes for injured vs healthy

# Create a heatmap for the top 100 genes
# First, create a new AnnData object with only the top 100 genes
top_genes_mask = adata_combined_subset.var_names.isin(top_genes_injured[:20])
top_genes_adata = adata_combined_subset[:, top_genes_mask]

# You may want to scale the data again just for this heatmap
sc.pp.scale(top_genes_adata)

# Plot the heatmap
sc.pl.heatmap(top_genes_adata, 
               var_names=top_genes_injured[:20], 
               groupby='condition', 
               cmap='viridis', 
               show_gene_labels=True)
plt.savefig(os.path.join(joined_output_base_dir, f"Heatmap_DEGs.png"), bbox_inches='tight')
plt.close()
# Convert DEGs to a list
deg_list = list(top_genes_injured)

# Perform GO term enrichment analysis
enr_go = gp.enrichr(gene_list=deg_list, 
                 gene_sets='GO_Biological_Process_2021',  # You can choose other gene sets too
                 organism='Mouse',  # or 'Human', depending on your data
                 outdir='test/enrichr_go')

# Show the results
print(enr_go.results.head(10))

go_results = enr_go.results

# Create a figure and axis
plt.figure(figsize=(10, 6))
plt.barh(go_results['Term'][:10], go_results['Adjusted P-value'][:10], color='skyblue')
plt.xlabel('Adjusted P-value')
plt.title('Top 10 GO Terms Enrichment')
plt.gca().invert_yaxis()  # Invert y-axis to have the highest p-value on top
plt.grid(axis='x')

# Save the plot
plt.savefig(os.path.join(joined_output_base_dir, 'go_enrichment_top_10_terms.png'), bbox_inches='tight')
plt.close()

# For KEGG pathway enrichment
enr_kegg = gp.enrichr(gene_list=deg_list, 
                 gene_sets='KEGG_2019_Mouse', 
                 organism='Mouse',
                 outdir='test/enrichr_kegg')

# Show the top pathways
print(enr_kegg.results.head(10))
kegg_results = enr_kegg.results

# Create a figure and axis
plt.figure(figsize=(10, 6))
plt.barh(kegg_results['Term'][:10], kegg_results['Adjusted P-value'][:10], color='lightgreen')
plt.xlabel('Adjusted P-value')
plt.title('Top 10 KEGG Pathways Enrichment')
plt.gca().invert_yaxis()  # Invert y-axis to have the highest p-value on top
plt.grid(axis='x')

# Save the plot
plt.savefig(os.path.join(joined_output_base_dir, 'kegg_enrichment_top_10_pathways.png'), bbox_inches='tight')
plt.close()  # Close the plot to avoid display if running in a notebook
'''
