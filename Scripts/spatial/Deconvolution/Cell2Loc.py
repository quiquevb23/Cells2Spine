'''
	We will do deconvolution of St data with T.Paralytica split into 2 conditions
'''
import os
import pandas as pd
from sklearn import metrics
import multiprocessing as mp
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import scanpy as sc
import scipy.sparse
import numpy as np
import cell2location

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

joined_base_dir = os.path.join(parent_dir, "Deconvolution")
#ref_signatures = os.path.join(joined_base_dir, "reference_signatures")
joined_output_base_dir = os.path.join(parent_output_dir, "Deconvolution")
run_name = f'{joined_base_dir}/cell2location_map'

os.makedirs(joined_base_dir, exist_ok=True)
#os.makedirs(ref_signatures, exist_ok=True)
os.makedirs(joined_output_base_dir, exist_ok=True)

single_cell_ref_h5 = "/storage/gge/Quique/TabulaeParalytica/single/GSE234774.h5"

adata_dict = {}


def deconvolution(adata_vis, sample, condition):
    mt_genes = adata_vis.var[adata_vis.var['mt'] == True].index
    adata_vis = adata_vis[:, ~adata_vis.var.index.isin(mt_genes)] # eliminate mt genes

    adata_vis.var_names_make_unique()
    # load regression model for sc for specific condition
    ref_signatures_path = os.path.join(joined_base_dir, f"reference_signatures_{condition}")
    adata_file = f"{ref_signatures_path}/sc.h5ad"
    adata_ref = sc.read_h5ad(adata_file)
    mod = cell2location.models.RegressionModel.load(f"{ref_signatures_path}", adata_ref)

    adata_ref = mod.export_posterior(
        adata_ref, use_quantiles=True,
        # choose quantiles
        add_to_varm=["q05","q50", "q95", "q0001"],
        #sample_kwargs={'batch_size': 2500, 'use_gpu': False}
    )

    # export estimated expression in each cluster
    if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
        inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_ref.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
    inf_aver.columns = adata_ref.uns['mod']['factor_names']
    print(inf_aver.iloc[0:5, 0:5])
    
    # find shared genes and subset both anndata and reference signatures
    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    # prepare anndata for cell2location model
    cell2location.models.Cell2location.setup_anndata(adata=adata_vis)
    # assign n of cells per location to 9 if healthy and 10 if injured
    l = 9 if condition == "healthy" else 10 if condition == "injured" else l
    # create and train the model
    mod = cell2location.models.Cell2location(
        adata_vis, cell_state_df=inf_aver,
        # the expected average cell abundance: tissue-dependent
        # hyper-prior which can be estimated from paired histology:
        N_cells_per_location=l,
        # hyperparameter controlling normalisation of
        # within-experiment variation in RNA detection:
        detection_alpha=20
    )
    mod.view_anndata_setup()
  
    mod.train(max_epochs=10000, 
              # train using full data (batch_size=None)
              batch_size=None,
              # use all data points in training because
              # we need to estimate cell abundance at all locations
              train_size=1,
              use_gpu=False,
              )

    # plot ELBO loss history during training, removing first 100 epochs from the plot
    mod.plot_history(1000)
    plt.legend(labels=['full data training']);
    plt.savefig(os.path.join(joined_output_base_dir, "ELBO_loss.png"))

    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata_vis = mod.export_posterior(
        adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': False}
    )

    # Save model
    mod.save(f"{run_name}/model_{sample}", overwrite=True)

    # mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

    # Save anndata object with results
    adata_file = f"{run_name}/sp{sample}.h5ad"
    adata_vis.write(adata_file)
    adata_file

    '''
    adata_file = f"{run_name}/sp.h5ad"
    adata_vis = sc.read_h5ad(adata_file)
    mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)
    '''
    
    adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
    #mod.plot_QC()

    # Loop through each unique label in adata.obs['_scvi_labels']
    for label in adata_vis.uns['mod']['factor_names']:
        #Generate spatial plot for the current cell type (label)
        sc.pl.spatial(adata_vis, 
                      cmap='magma',          # Set the color map
                      color=label,           # Use the cell type (label) as the color
                      ncols=5,               # Number of columns for subplot layout
                      size=1.3,              # Adjust the size of the spots
                      img_key='hires',       # Use the high-resolution image
                      vmin=0,                # Minimum color value
                      vmax='p99.2',          # Maximum color value at the 99.2 percentile
                      show=False             # Set to False if you don’t want to display each plot immediately
                     )

        # Save each plot with a distinct name using the current label
        plt.savefig(os.path.join(joined_output_base_dir, f"deconvolution_{sample}_{label}.png"), dpi=300, bbox_inches='tight')
        plt.close()

    from cell2location import run_colocation
    res_dict, adata_vis = run_colocation(
        adata_vis,
        model_name='CoLocatedGroupsSklearnNMF',
        train_args={
          'n_fact': np.arange(11, 13), # IMPORTANT: use a wider range of the number of factors (5-30)
          'sample_name_col': 'sample', # columns in adata_vis.obs that identifies sample
          'n_restarts': 3 # number of training restarts
        },
        # the hyperparameters of NMF can be also adjusted:
        model_kwargs={'alpha': 0.01, 'init': 'random', "nmf_kwd_args": {"tol": 0.000001}},
        export_args={'path': f'{run_name}/CoLocatedComb_{sample}/'}
    )
    # Plot NMF weights to now co-location of cell types
    res_dict['n_fact12']['mod'].plot_cell_type_loadings()
    plt.savefig(os.path.join(joined_output_base_dir, f"heatmap_colocation_{sample}.png"), dpi=300, bbox_inches='tight')
    # Compute expected expression per cell type
    expected_dict = mod.module.model.compute_expected_per_cell_type(
        mod.samples["post_sample_q05"], mod.adata_manager
    )
    
    # Add to anndata blayers
    for i, n in enumerate(mod.factor_names_):
        adata_vis.layers[n] = expected_dict['mu'][i]

    # Save anndata object with results
    adata_file = f"{run_name}/sp.h5ad"
    adata_vis.write(adata_file)

def create_reference_model(adata_ref, condition):
    #we will create ref signature model for each condition
    ref_signatures = os.path.join(joined_base_dir, f"reference_signatures_{condition}")
    os.makedirs(ref_signatures, exist_ok=True)
    # Rename duplicates
    adata_ref.var_names_make_unique()
    adata_ref.obs_names_make_unique()

    del adata_ref.raw
    from cell2location.utils.filtering import filter_genes
    selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

    # filter the object
    adata_ref = adata_ref[:, selected].copy()
    # prepare anndata for the regression model
    cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                            # batch they come from
                            batch_key='label',
                            # cell type, covariate used for constructing signatures
                            labels_key='cell_l4')
                            # multiplicative technical effects (platform, 3' vs 5', donor effect)
                            # categorical_covariate_keys=['Method']

    # create the regression model
    from cell2location.models import RegressionModel
    mod = RegressionModel(adata_ref)

    # view anndata_setup as a sanity check
    mod.view_anndata_setup()
    mod.train(max_epochs=250, use_gpu=False)
    #mod.plot_history(20)
    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata_ref = mod.export_posterior(
        adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': False}
    )

    # Save model
    mod.save(ref_signatures, overwrite=True)

    # Save anndata object with results
    adata_file = os.path.join(ref_signatures, "sc.h5ad")
    adata_ref.write(adata_file)
    print(adata_file)


# First create the regression mod for reference single cell
##############################

adata_ref = sc.read(single_cell_ref_h5)
# Split adata into 2 based on injured (7d, 14d) and healthy individuals:

adata_ref_healthy = adata_ref[adata_ref.obs['label'].isin(['uninjured'])].copy()
adata_ref_injured = adata_ref[adata_ref.obs['label'].isin(['7d', '14d'])].copy()

condition_healthy = "healthy"
condition_injured = "injured"

create_reference_model(adata_ref_healthy, condition_healthy)
create_reference_model(adata_ref_injured, condition_injured)

'''
#Load the model (done in function)
adata_file = os.path.join(ref_signatures,"sc.h5ad")
adata_ref = sc.read_h5ad(adata_file)
mod = cell2location.models.RegressionModel.load(f"{ref_signatures}", adata_ref)


adata_ref = mod.export_posterior(
    adata_ref, use_quantiles=True,
    # choose quantiles
    add_to_varm=["q05","q50", "q95", "q0001"],
#    sample_kwargs={'batch_size': 2500, 'use_gpu': False}
)
#mod.plot_QC()
'''

########################################################

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
#            adata = st.convert_scanpy(adata)
#            adata_dict[sample_name] = stSME_normalization(adata, sample_name)
            condition = adata.obs['condition'].unique()[0]
            deconvolution(adata, sample_name, condition)

