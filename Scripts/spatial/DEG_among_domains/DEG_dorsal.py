"""
	Subset dorsal region of adata for DGE among healthy and injured
"""

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

base_dir = args.base_dir
output_base_dir = args.output_base_dir

#Create new folders for "joined" datasets by Harmony
parent_dir = os.path.dirname(base_dir)
parent_output_dir = os.path.dirname(output_base_dir)

joined_base_dir = os.path.join(parent_dir, "DEGs")
joined_output_base_dir = os.path.join(parent_output_dir, "DEGs")

os.makedirs(joined_base_dir, exist_ok=True)
os.makedirs(joined_output_base_dir, exist_ok=True)

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
    #sc.pp.scale(adata)
    return adata

def get_celltype_proportions(adata):
    cell2loc_file_path = "/storage/gge/Quique/Cells2SpineData/Pilot/spatial/matrices/Deconvolution/cell2location_map"
    adata_cell2loc_3 = sc.read(os.path.join(cell2loc_file_path, "spSpatial_3.h5ad"))
    adata_cell2loc_4 = sc.read(os.path.join(cell2loc_file_path, "spSpatial_4.h5ad"))

    adata_cell2loc_3.obs[adata_cell2loc_3.uns['mod']['factor_names']] = adata_cell2loc_3.obsm['q05_cell_abundance_w_sf']
    adata_cell2loc_4.obs[adata_cell2loc_4.uns['mod']['factor_names']] = adata_cell2loc_4.obsm['q05_cell_abundance_w_sf']

    # Subset for the specified cluster in both datasets
    subset_3 = adata_cell2loc_3[adata_cell2loc_3.obs['leiden_0.4'] == '6']
    subset_4 = adata_cell2loc_4[adata_cell2loc_4.obs['leiden_0.4'] == '6']
    
    # Calculate mean proportions for each cell type
    mean_proportions_3 = subset_3.obs[subset_3.uns['mod']['factor_names']].mean()
    mean_proportions_4 = subset_4.obs[subset_4.uns['mod']['factor_names']].mean()

    mean_of_means = (mean_proportions_3 + mean_proportions_4) / 2
    proportions_df = pd.DataFrame(mean_of_means, columns=['Mean_Proportion']).reset_index()
    proportions_df.columns = ['Cell_Type', 'Mean_Proportion']
    proportions_df.to_csv(os.path.join(joined_output_base_dir, 'proportions_dorsal_injured.csv'))

def process(adata_combined):
    condition = adata_combined.obs['condition'].unique()[0]
    if condition == 'injured':
        subset_adata = adata_combined[adata_combined.obs['leiden_0.4'] == '6'].copy() #leiden cluste$
        sc.pp.normalize_total(subset_adata, target_sum=1e4)
        
        #subset_adata = normalize(subset_adata)
        return subset_adata
    else:
        subset_adata = adata_combined[adata_combined.obs['highlight_group'] == 'highlight'].copy() #copy selection for dorsal site of healthy
        sc.pp.normalize_total(subset_adata, target_sum=1e4)

        #subset_adata = normalize(subset_adata)
        return subset_adata

#def biomart_convert():


def find_unique_genes(healthy_subset, injured_subset, threshold=0.1):
    """
    Identify genes uniquely expressed in either the healthy or injured condition.

    
    Parameters:
    -----------
    healthy_subset : AnnData
        Subset of dorsal region for the healthy condition.
    injured_subset : AnnData
        Subset of dorsal cluster for the injured condition.
    threshold : float
        Expression threshold (in log-transformed data) to determine whether a gene is considered "present".
        Typically, a threshold of around 0.1 - 1 is appropriate for log-transformed data.
    
    """
    # Filter out genes not present in at least 3 spots to avoid ties in log2fold change plot
    sc.pp.filter_genes(healthy_subset, min_cells=3)
    sc.pp.filter_genes(injured_subset, min_cells=3)

    # Get all unique gene names from both healthy and injured datasets
    healthy_genes = set(healthy_subset.var_names)
    injured_genes = set(injured_subset.var_names)
    all_genes = list(healthy_genes.union(injured_genes))
    
    # Initialize dictionaries to store mean expression values
    healthy_expression = {}
    injured_expression = {}
    
    # Fill in mean expression for healthy condition (set missing genes to 0)
    for gene in healthy_genes:
        healthy_expression[gene] = healthy_subset[:, gene].X.mean()
    
    # Fill in mean expression for injured condition (set missing genes to 0)
    for gene in injured_genes:
        injured_expression[gene] = injured_subset[:, gene].X.mean()
    # Convert to pandas Series for easy comparison
    healthy_expression = pd.Series(healthy_expression)
    injured_expression = pd.Series(injured_expression)
    
    # Find genes uniquely expressed in the healthy condition
    unique_healthy_genes = healthy_expression[
        (healthy_expression > threshold) & (injured_expression <= threshold)
    ].index.tolist()
    
    # Find genes uniquely expressed in the injured condition
    unique_injured_genes = injured_expression[
        (injured_expression > threshold) & (healthy_expression <= threshold)
    ].index.tolist()

    ##Now let's plot the unique genes for each condition
    healthy_counts = healthy_subset[:, unique_healthy_genes].X
    injured_counts = injured_subset[:, unique_injured_genes].X
    # Convert counts to DataFrame
    healthy_df = pd.DataFrame(healthy_counts, columns=unique_healthy_genes, index=healthy_subset.obs.index)
    injured_df = pd.DataFrame(injured_counts, columns=unique_injured_genes, index=injured_subset.obs.index)
    # Calculate mean counts for unique genes
    healthy_means = healthy_df.mean().reset_index()
    injured_means = injured_df.mean().reset_index()
    # Prepare the DataFrame for CSV export
    combined_unique_genes_df = pd.concat([
        healthy_means.assign(condition='Healthy'),
        injured_means.assign(condition='Injured')
    ], ignore_index=True)
    
    combined_unique_genes_df.columns = ['Gene', 'Mean Count', 'Condition']
    combined_unique_genes_df.to_csv(os.path.join(joined_output_base_dir, 'unique_genes_counts.csv'), index=False)


    # Step 3: Prepare the data for plotting
    plt.figure(figsize=(12, 6))
    
    # Box plot for unique genes
    plt.subplot(1, 2, 1)
    sns.boxplot(data=healthy_df, orient='h')
    plt.title('Unique Genes in Healthy Condition')
    plt.xlabel('Counts')
    plt.yticks([])  # Remove y-axis labels

    plt.subplot(1, 2, 2)
    sns.boxplot(data=injured_df, orient='h')
    plt.title('Unique Genes in Injured Condition')
    plt.xlabel('Counts')
    plt.yticks([])  # Remove y-axis labels

    plt.tight_layout()
    plt.savefig(os.path.join(joined_output_base_dir,'unique_genes_counts_plot.png'), bbox_inches='tight')  # Save plot as PNG 
    
    return unique_healthy_genes, unique_injured_genes

def plot_logfold_change(deg_results):
    """
    Function to plot the logFold Change of DEGs among conditions
    """
    pval_threshold = 0.05  # Adjusted p-value threshold
    logfc_threshold = 1.0   # Log fold change threshold

    # Identify top 1000 most significant genes by p-value
    top_genes_by_pval = deg_results.nsmallest(1000, 'pvals_adj')

    # Filter for significant DEGs
    filtered_degs = deg_results[
        (deg_results['pvals_adj'] < pval_threshold) &
        (deg_results['logfoldchanges'].abs() > logfc_threshold)
    ]

    # Create a volcano plot for visualization
    plt.figure(figsize=(10, 6))
    sns.scatterplot(
        data=deg_results,
        x='logfoldchanges',
        y=-np.log10(deg_results['pvals_adj']),
        alpha=0.7,
        color='grey'
    )

    # Highlight significant DEGs
    sns.scatterplot(
        data=filtered_degs,
        x='logfoldchanges',
        y=-np.log10(filtered_degs['pvals_adj']),
        color='red'
    )

    # Highlight top 1000 significant genes in orange
    sns.scatterplot(
        data=top_genes_by_pval,
        x='logfoldchanges',
        y=-np.log10(top_genes_by_pval['pvals_adj']),
        color='orange',
        label='Top 1000 by P-Value'
    )

    # Add labels and title
    plt.axhline(y=-np.log10(0.05), color='blue', linestyle='--')  # Significance threshold line
    plt.title('Volcano Plot of DEGs')
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 Adjusted P-Value')

    plt.autoscale()
    #plt.xlim(-5, 5)  # Adjust limits based on your data
    #plt.ylim(0, 10)  # Adjust limits based on your data

    plt.savefig(os.path.join(joined_output_base_dir, "volcano_plot_degs.png"), bbox_inches='tight')
    plt.close()


def DEG_dorsal(adata_combined, unique_healthy_genes, unique_injured_genes):
    print("Processing dorsal subset")

#    adata_combined_subset = normalize(adata_combined_subset)
    # Because each subset is already normalized and log-transformed, we do not need to renormalize
    unique_genes = set(unique_healthy_genes + unique_injured_genes)
    all_genes = set(adata_combined.var_names)
    non_unique_genes = list(all_genes.difference(unique_genes))

    # Subset to include only non-unique genes
    adata_combined_subset = adata_combined[:, non_unique_genes].copy()

    adata_combined_subset = normalize(adata_combined_subset)
    # Perform differential expression analysis between the two groups: if not normality better use wilcoxon
    sc.tl.rank_genes_groups(adata_combined_subset, groupby='condition', groups=['injured'], reference='healthy', method='wilcoxon', use_raw=False)

    # View the results
    sc.pl.rank_genes_groups(adata_combined_subset, n_genes=20, sharey=False)
    plt.savefig(os.path.join(joined_output_base_dir, f"Rank_genes_groups.png"), bbox_inches='tight')
    plt.close()

    # Extract the results for the injured group
    deg_results = adata_combined_subset.uns['rank_genes_groups']

    # Convert to pandas DataFrame for the injured group
    deg_df = pd.DataFrame({
        'gene': deg_results['names']['injured'],
        'logfoldchanges': deg_results['logfoldchanges']['injured'],
        'pvals': deg_results['pvals']['injured'],
        'pvals_adj': deg_results['pvals_adj']['injured'],
        'scores': deg_results['scores']['injured']
    })

    plot_logfold_change(deg_df)

    # Save the top 100 DE genes for injured group to a CSV file
    csv_output_path = os.path.join(joined_output_base_dir, "Top_DEG_genes_injured_vs_healthy.csv")
    deg_df.to_csv(csv_output_path, index=False)

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



    #Before doing GO term, need to convert to mouse
    '''
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

    # The following is for doing only ORA on samples without DGE
    '''
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
    '''
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

datas_healthy = []
datas_injured = []
adata_dict = {}
datas_dorsal = []

# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name, "outs", "matrices")
    for file in os.listdir(sample_path):
        if file.endswith('_highlight.h5ad'):
            sample_name = file_name
            print(f"Processing healthy adata {sample_name}")
            # Define the directories to save plots and matrices
            output_dir = os.path.join(output_base_dir, sample_name)
            os.makedirs(output_dir, exist_ok=True)
            # This has the annotation for clusters but is not normalized

            adata = sc.read_h5ad(os.path.join(sample_path, file))
            #adata = st.convert_scanpy(adata)
            #adata_dict[sample_name] = stSME_normalization(adata, sample_name)
            adata_dict[sample_name] = adata
            datas_healthy.append(adata_dict[sample_name])

        elif file.startswith('adata_stlearn_common_domains_'):
            sample_name = file_name
            # Define the directories to save plots and matrices
            output_dir = os.path.join(output_base_dir, sample_name)
            os.makedirs(output_dir, exist_ok=True)
            # This has the annotation for clusters but is not normalized

            adata = sc.read_h5ad(os.path.join(sample_path, file))
            condition = adata.obs['condition'].unique()[0]
            if condition == 'injured':
                #adata = st.convert_scanpy(adata)
                print(f"Processing injured adata {sample_name}")
                #adata_dict[sample_name] = stSME_normalization(adata, sample_name)
                adata_dict[sample_name] = adata
                datas_injured.append(adata_dict[sample_name])


adata_healthy = sc.concat(datas_healthy, join='inner', index_unique = "_")
adata_injured = sc.concat(datas_injured, join='inner', index_unique = "_")
#Process to get only dorsal regions
healthy_subset = process(adata_healthy)
injured_subset = process(adata_injured)

datas_dorsal.append(healthy_subset)
datas_dorsal.append(injured_subset)

dorsal = sc.concat(datas_dorsal, join='inner', index_unique = "_")

#Find unique genes in each condition
'''
unique_healthy_genes, unique_injured_genes = find_unique_genes(healthy_subset, injured_subset, threshold=0.1) 

datas_dorsal.append(healthy_subset)
datas_dorsal.append(injured_subset)

dorsal = sc.concat(datas_dorsal, join='inner', index_unique = "_")
#Perform DEA between dorsal region
DEG_dorsal(dorsal, unique_healthy_genes, unique_injured_genes)
'''
