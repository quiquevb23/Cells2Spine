'''
	Script to perform individual and joint clustering with StLearn
'''
import os
import pandas as pd
from sklearn import metrics
import multiprocessing as mp
import matplotlib.pyplot as plt
import argparse
import stlearn as st
import scanpy as sc
import scipy.sparse
import numpy as np
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

joined_base_dir = os.path.join(base_dir, "joined_cond_all", "stlearn_dir")
joined_output_base_dir = os.path.join(output_base_dir, "joined_cond_all", "stlearn_dir")

os.makedirs(joined_base_dir, exist_ok=True)
os.makedirs(joined_output_base_dir, exist_ok=True)

def stlearn_clustering(data, sample, output_dir, sample_path):
    tile_path = os.path.join("./tmp/tiles", sample)
    os.makedirs(tile_path, exist_ok=True)

    # Create a layer named stlearn to store changes
    if 'stlearn' not in data.layers:
        data.layers['stlearn'] = data.X.copy()
        print("Created 'stlearn' layer from adata.X.")

    # Switch adata.X to the 'stlearn' layer for processing 
    data.X = data.layers['stlearn']
    # Processing of counts data
    st.pp.filter_genes(data,min_cells=1)
    data.raw = data.copy() # copy here to have the same nre of variables after filtering lowly expressed genes
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

    output_file = os.path.join(sample_path, f"{sample}_stlearn.h5ad")
    data_SME.write(output_file) #save the stlearn files 

#    return data_SME #return data here without scaling if we want to norm on aggr dataset (approach3)

    st.pp.scale(data_SME)
    st.em.run_pca(data_SME,n_comps=50)
    # K-means clustering on stSME normalised PCA
    st.tl.clustering.kmeans(data_SME,n_clusters=12, use_data="X_pca", key_added="X_pca_kmeans")
    st.pl.cluster_plot(data_SME, use_label="X_pca_kmeans", size=20)
    plt.savefig(os.path.join(output_dir, f"clustering_kmeans_{sample}.png"))
    # louvain clustering on stSME normalised data
    st.pp.neighbors(data_SME,n_neighbors=17,use_rep='X_pca')
    st.tl.clustering.louvain(data_SME, resolution=1.19)
    st.pl.cluster_plot(data_SME,use_label="louvain", size=20)
    plt.savefig(os.path.join(output_dir, f"clustering_louvain_{sample}.png"))
    
    
    return data_SME #return data after scaling if we are not going to normalize on aggr (approach2)


def process_combined(adata_combined):
    # Normalize concatenated object 
    condition = adata_combined.obs['condition'].unique()[0]
    sc.pp.filter_genes(adata_combined, min_cells=3) #filter genes
#    sc.pp.filter_cells(adata_combined, min_counts=1) #filter genes

#### Comment this block to follow approach 2 / uncomment for approach 3
    '''
    # Normalize data
    sc.pp.normalize_total(adata_combined, target_sum=1e4)
    # Log transformation
    sc.pp.log1p(adata_combined)

    # Store raw data
    adata_combined.raw = adata_combined

    # Calculate HVGs
    sc.pp.highly_variable_genes(adata_combined, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata_combined = adata_combined[:, adata_combined.var['highly_variable']].copy()

    sc.pp.scale(adata_combined, max_value=10) # only scale
    '''
####
    sc.pp.pca(adata_combined, n_comps=30, svd_solver='arpack')
    sc.pl.pca_scatter(adata_combined, color=["sample", "condition"])
    plt.savefig(os.path.join(joined_output_base_dir, f"PCA_scatter.png"), bbox_inches='tight')
    plt.close()

    sc.pl.pca_variance_ratio(adata_combined)
    plt.savefig(os.path.join(joined_output_base_dir, f"PCA_variance_ratio.png"), bbox_inches='tight')
    plt.close()

    '''
    # We will compute now UMAP on combined dataset without integration
    adata_noint = adata_combined.copy()
    sc.pp.neighbors(adata_noint)
    sc.tl.leiden(adata_noint, resolution=0.4, key_added='leiden_0.4')
    sc.tl.umap(adata_noint)
    sc.pl.umap(adata_noint, color=["sample", "condition", "leiden_0.4"], wspace=0.5)
    plt.savefig(os.path.join(joined_output_base_dir, f"UMAP_unintegrated.png"), bbox_inches='tight')
    plt.close()
    adata_noint.write(os.path.join(joined_base_dir, f"adata_no_integrated.h5ad"))
    '''

    # Now run Harmony
    # Since we have normalized and scaled data after concat, we can run now Harmony, with first 50 P$

    meta_data = adata_combined.obs
    data_mat = adata_combined.obsm['X_pca']
    import harmonypy as hm
    ho = hm.run_harmony(data_mat, meta_data, 'sample')
    adata_combined.obsm['X_pca'] = ho.Z_corr.T

    sc.pp.neighbors(adata_combined, n_pcs=30)
    sc.tl.umap(adata_combined)
    sc.tl.leiden(adata_combined, resolution=0.4, key_added='leiden_0.4')
    sc.pl.umap(adata_combined, color=["sample", "condition", "leiden_0.4"], wspace=0.5)
    plt.savefig(os.path.join(joined_output_base_dir, f"UMAP_integrated_{condition}.png"), bbox_inches='tight')
    plt.close()
    adata_combined.write(os.path.join(joined_base_dir, f"adata_integrated.h5ad"))

    return adata_combined

def process_all(adata_combined):
    # Normalize concatenated object 
#### Comment this block to follow approach 2 / uncomment for approach 3
    '''
    # Normalize data
    sc.pp.normalize_total(adata_combined, target_sum=1e4)
    # Log transformation
    sc.pp.log1p(adata_combined)

    # Store raw data
    adata_combined.raw = adata_combined

    # Calculate HVGs
    sc.pp.highly_variable_genes(adata_combined, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata_combined = adata_combined[:, adata_combined.var['highly_variable']].copy()

    sc.pp.scale(adata_combined, max_value=10) # only scale
    '''
####
    #Then we plot the clustes spatially on each slice
    sc.pp.pca(adata_combined, n_comps=30, svd_solver='arpack')
    # Now run Harmony
    # Since we have normalized and scaled data after concat, we can run now Harmony, with first 50 P$

    meta_data = adata_combined.obs
    data_mat = adata_combined.obsm['X_pca']
    import harmonypy as hm
    ho = hm.run_harmony(data_mat, meta_data, 'condition')
    adata_combined.obsm['X_pca'] = ho.Z_corr.T

    sc.pp.neighbors(adata_combined, n_pcs=30)
    sc.tl.umap(adata_combined)
    sc.tl.leiden(adata_combined, resolution=0.4, key_added='leiden_0.4')
    sc.tl.leiden(adata_combined, resolution=0.8, key_added='leiden_0.8')
    # Create color map paletter
    sample_unique = adata_combined.obs['sample'].unique().tolist()
    condition_unique = adata_combined.obs['condition'].unique().tolist()
    leiden_0_4_unique = adata_combined.obs['leiden_0.4'].cat.categories.tolist()
    leiden_0_8_unique = adata_combined.obs['leiden_0.8'].cat.categories.tolist()

    # Step 2: Create distinct palettes for leiden_0.4 and leiden_0.8
    # Use different Seaborn color palettes or any color generator
    leiden_0_4_palette = sns.color_palette("bright", len(leiden_0_4_unique))
    leiden_0_8_palette = sns.color_palette("Paired", len(leiden_0_8_unique))
    sample_palette = sns.color_palette("tab10", len(sample_unique))
    condition_palette = sns.color_palette("Dark2", len(condition_unique))

    # We zip the cluster labels with the generated colors
    leiden_0_4_color_mapping = {category: leiden_0_4_palette[i] for i, category in enumerate(leiden_0_4_unique)}
    leiden_0_8_color_mapping = {category: leiden_0_8_palette[i] for i, category in enumerate(leiden_0_8_unique)}
    sample_color_mapping = {category: sample_palette[i] for i, category in enumerate(sample_unique)}
    condition_color_mapping = {category: condition_palette[i] for i, category in enumerate(condition_unique)}

    combined_palette = {
        **sample_color_mapping,
        **condition_color_mapping,
        **leiden_0_4_color_mapping,
        **leiden_0_8_color_mapping
    }

    sc.pl.umap(adata_combined, color="sample", palette=sample_color_mapping, wspace=0.5, show=False)
    plt.savefig(os.path.join(joined_output_base_dir, f"UMAP_sample.png"), bbox_inches='tight')
    plt.close()

    sc.pl.umap(adata_combined, color="condition", palette=condition_color_mapping, wspace=0.5, show=False)
    plt.savefig(os.path.join(joined_output_base_dir, f"UMAP_condition.png"), bbox_inches='tight')
    plt.close()

    sc.pl.umap(adata_combined, color="leiden_0.4", palette=leiden_0_4_color_mapping, wspace=0.5, show=False)
    plt.savefig(os.path.join(joined_output_base_dir, f"UMAP_leiden_0_4.png"), bbox_inches='tight')
    plt.close()

    sc.pl.umap(adata_combined, color="leiden_0.8", palette=leiden_0_8_color_mapping, wspace=0.5, show=False)
    plt.savefig(os.path.join(joined_output_base_dir, f"UMAP_leiden_0_8.png"), bbox_inches='tight')
    plt.close()

    adata_combined.write(os.path.join(joined_base_dir, f"adata_integrated.h5ad"))

    for file_name in os.listdir(os.path.join(base_dir, "indiv_samples")):
        sample_path = os.path.join(base_dir, "indiv_samples", file_name, "outs", "matrices")
        for file in os.listdir(sample_path):
            if file.endswith('stlearn.h5ad'):
                sample_name = file_name
                
                # Define the directories to save plots and matrices
                output_dir = os.path.join(output_base_dir, sample_name)
                os.makedirs(output_dir, exist_ok=True)

                adata = sc.read_h5ad(os.path.join(sample_path, file))

                filtered_adata = adata_combined[adata_combined.obs['sample'] == sample_name].copy()

                adata.obs['leiden_0.4'] = pd.Categorical(filtered_adata.obs['leiden_0.4'].values)
                adata.obs['leiden_0.8'] = pd.Categorical(filtered_adata.obs['leiden_0.8'].values)

                # Set the color palette for the 'domain' column
                sc.pl.spatial(adata, color=['leiden_0.4'], wspace=0.5, palette=leiden_0_4_color_mapping)
                plt.savefig(os.path.join(joined_output_base_dir, f'leiden_clusters_{sample_name}.png'),  bbox_inches='tight')
                plt.close()
                
                # This adata will always contain the attribute .raw so if want we can:
                # adata.X = adata.raw.X.copy() to get back to raw counts
                adata.write(os.path.join(sample_path, f"adata_stlearn_common_domains_{sample_name}.h5ad"))

datas_healthy = []
datas_injured = []

adata_dict = {}
for file_name in os.listdir(os.path.join(base_dir, "indiv_samples")):
    sample_path = os.path.join(base_dir, "indiv_samples", file_name, "outs", "matrices")
    for file in os.listdir(sample_path):
        if file.endswith('qc_metrics.h5ad'): # we take the raw data
            sample_name = file_name

            # Define the directories to save plots and matrices
            output_dir = os.path.join(output_base_dir, sample_name, "stlearn_dir")
            os.makedirs(output_dir, exist_ok=True)

            adata = sc.read_h5ad(os.path.join(sample_path, file))
            adata = st.convert_scanpy(adata)
            print(f"Processing sample {sample_name} with size: {adata.shape}")
            adata.obs['sample'] = file_name
            # Create raw attribute set it to adata.X (raw and unnormalized)
            if adata.raw is None:
                adata.raw = adata.copy()
                print("Initialized adata.raw with a copy of the current AnnData.")
            #If we want to follow approach 1, we have to concatenate the 4 samples without normalizing
            if file_name == "Spatial_1" or file_name == "Spatial_2":
                adata.obs["condition"] = "healthy"
                adata_dict[sample_name] = stlearn_clustering(adata, sample_name, output_dir, sample_path) # perform clustering and return adata
                #adata_dict[sample_name].write(os.path.join(sample_path, f"adata_stlearn_domains_{sample_name}.h5ad"))
                datas_healthy.append(adata_dict[sample_name]) #add to healthy adatas after SME norma$
                # datas.append(adata) #append without processing raw (approach 1)
            else:
                adata.obs["condition"] = "injured"
                adata_dict[sample_name] = stlearn_clustering(adata, sample_name, output_dir, sample_path) # perform clustering and return adata
                #adata_dict[sample_name].write(os.path.join(sample_path, f"adata_stlearn_domains_{sample_name}.h5ad"))
                datas_injured.append(adata_dict[sample_name])

                # datas.append(adata) #append without processing raw (approach 1) 

#concatenate only with common genes using join inner
'''
	Depends on what we want to do: if we join by outer, the aggregated datasets can have 0s (if gene was not present originally). In this case we should use
	a method like Scran that can handle 0s. 
	If we want to use normalize_total after concatenating, we should concatenate with inner (default), to ensure only common genes are present in dataset
	In any case, we should not scale individual datasets since this may make HVG to not work
'''

adata_healthy = sc.concat(datas_healthy, join='inner', index_unique = "_")
adata_injured = sc.concat(datas_injured, join='inner', index_unique = "_")

#Process per condition with Harmony first
adata_healthy = process_combined(adata_healthy)
adata_injured = process_combined(adata_injured)

datas_all = []
datas_all.append(adata_healthy)
datas_all.append(adata_injured)

adata_all = sc.concat(datas_all, join='inner', index_unique = "_")

#Process combined
process_all(adata_all)

