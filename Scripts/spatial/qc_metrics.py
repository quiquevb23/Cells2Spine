"""
	Script to check QC metrics plots and decide which filters to apply for removing low-quality
        spots - specifically first those that fall out of the tissue
"""
import scanpy as sc
import anndata as ad
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import median_abs_deviation
import argparse
import squidpy

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

# We will identify cells as outliers if they differ by 5 MADs for the QC in particular

def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name, "outs")
    for file in os.listdir(sample_path):
        if file.endswith('filtered_feature_bc_matrix.h5'):
            sample_name = file_name
            print(f"Processing sample {sample_name}")
            # Define the directories to save plots and matrices
            output_dir = os.path.join(output_base_dir, sample_name)
            output_data_dir = os.path.join(base_dir, file_name, "outs", "matrices")
            os.makedirs(output_dir, exist_ok=True)
            os.makedirs(output_data_dir, exist_ok=True)
            #Read Visium matrices
            counts_file = os.path.join(sample_path, file)
            adata = squidpy.read.visium(path=sample_path, counts_file=counts_file) #specify path to Visium data
            adata.var_names_make_unique()
            adata.layers['counts'] = adata.X #save raw counts
            adata.var['mt'] = adata.var_names.str.startswith('mt-')
            adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
            adata.var["hb"] = adata.var_names.str.startswith(("Hbb", "Hba", "Hbq"))
            sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo','hb'],
                log1p=True, inplace=True)
            # QC metrics are calculated with unnormalized data and stored, so if we want to use them later they are unnormalized

            sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],
              jitter=0.4, multi_panel=True, show=False)
            plt.savefig(os.path.join(output_dir, 'violin.png'), bbox_inches='tight')

            plt.close()
            # Plot scatter and distplots
            sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color='pct_counts_mt', show=False)
            plt.savefig(os.path.join(output_dir, 'scatter.png'), bbox_inches='tight')
            plt.close()
            # Distplots
            sns.displot(adata.obs['total_counts'], bins=100, kde=False)
            plt.savefig(os.path.join(output_dir, 'total_counts_displot.png'), bbox_inches='tight')
            plt.close()
            sns.displot(adata.obs['pct_counts_mt'], bins=100, kde=False)
            plt.savefig(os.path.join(output_dir, 'pct_counts_mt_displot.png'), bbox_inches='tight')
            plt.close()
            # Extract summary statistics
            df = adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"]]
            summary_stats = df.describe().T  # Get summary statistics (mean, std, min, 25%, 50%, 75%$

            '''
            # Get additional specific statistics
            specific_stats = df.agg(['mean', 'median', 'quantile']).T
            specific_stats['quantile_25'] = df.quantile(0.25)
            specific_stats['quantile_75'] = df.quantile(0.75)
            '''

            # Plot scatter and distplot centered in median
            # Get median and IQR for filtering
            median_total_counts = df['total_counts'].median()
            median_n_genes = df['n_genes_by_counts'].median()
            iqr_total_counts = df['total_counts'].quantile(0.75) - df['total_counts'].quantile(0.25)
            iqr_n_genes = df['n_genes_by_counts'].quantile(0.75) - df['n_genes_by_counts'].quantile(0.25)
            
            # Filter the data within 1.5 IQR of the median
            filtered_data = df[
                (df['total_counts'] > (median_total_counts - 1.5 * iqr_total_counts)) &
                (df['total_counts'] < (median_total_counts + 1.5 * iqr_total_counts)) &
                (df['n_genes_by_counts'] > (median_n_genes - 1.5 * iqr_n_genes)) &
                (df['n_genes_by_counts'] < (median_n_genes + 1.5 * iqr_n_genes))
            ]
            # Scatter plot with filtered data
            adata_filtered = adata[filtered_data.index]

            # Create the scatter plot centered around the median
            sc.pl.scatter(adata_filtered, x='total_counts', y='n_genes_by_counts', color='pct_counts_mt', 
                          title='Filtered Scatter Plot (centered around median)', show=False)
            plt.savefig(os.path.join(output_dir, 'scatter_centered_median.png'), bbox_inches='tight')
            plt.close()

            # Distplot centered in median
            sns.displot(adata_filtered.obs['total_counts'], bins=100, kde=False)
            plt.savefig(os.path.join(output_dir, 'total_counts_centered_median_displot.png'), bbox_inches='tight')
            plt.close()
            sns.displot(adata_filtered.obs['pct_counts_mt'], bins=100, kde=False)
            plt.savefig(os.path.join(output_dir, 'pct_counts_mt_centered_median_displot.png'), bbox_inches='tight')
            plt.close()
            #PCT counts MT centered in median is taking into account median for total counts and n_genes by counts

            # We will now define if a cell is outlier
            adata.obs["outlier"] = (
                is_outlier(adata, "log1p_total_counts", 5)
                | is_outlier(adata, "log1p_n_genes_by_counts", 5)
            )
            # We will filter cells with MADs for pct_counts_mt more 3, and cells with > 8% MT genes
            adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3)
          #      adata.obs["pct_counts_mt"] > 8 #this can be excluded

            adata.obs["ribo_outlier"] = is_outlier(adata, "pct_counts_ribo", 3)

            # Define observations if they fall below Q25, Q50, and Q75 for 'total_counts'
            # Add a column to indicate total counts below the 25th percentile, median, and 75th percentile
            quantile_25 = df['total_counts'].quantile(0.25)
            median = df['total_counts'].median()
            quantile_75 = df['total_counts'].quantile(0.75)

            adata.obs['total_counts_category'] = pd.cut(adata.obs['total_counts'],
                                                       bins=[-float('inf'), quantile_25, median, quantile_75, float('inf')],
                                                       labels=['below_q25', 'below_median', 'below_q75', 'above_q75'])
            # Segregate high and low outliers for each metric
            def get_high_low_non_outliers(adata, metric, nmads):
                outliers = is_outlier(adata, metric, nmads)
                non_outliers = adata.obs[~outliers]
                high_value = non_outliers[metric].max()
                low_value = non_outliers[metric].min()
                return high_value, low_value

            # Get high and low values for non-outliers for each metric
            pct_counts_mt_high, pct_counts_mt_low = get_high_low_non_outliers(adata, "pct_counts_mt", 3)
            pct_counts_ribo_high, pct_counts_ribo_low = get_high_low_non_outliers(adata, "pct_counts_ribo", 3)
            log1p_total_counts_high, log1p_total_counts_low = get_high_low_non_outliers(adata, "log1p_total_counts", 5)
            log1p_n_genes_by_counts_high, log1p_n_genes_by_counts_low = get_high_low_non_outliers(adata, "log1p_n_genes_by_counts", 5)
            total_counts_high, total_counts_low = get_high_low_non_outliers(adata, "total_counts", 5)
            n_genes_by_counts_high, n_genes_by_counts_low = get_high_low_non_outliers(adata, "n_genes_by_counts", 5)

            # Get high and low values for outliers for each metric
            pct_counts_mt_high_outlier = adata.obs[adata.obs['mt_outlier']==True]['pct_counts_mt'].max()
            pct_counts_mt_low_outlier = adata.obs[adata.obs['mt_outlier']==True]['pct_counts_mt'].min()
            pct_counts_ribo_high_outlier = adata.obs[adata.obs['ribo_outlier']==True]['pct_counts_ribo'].max()
            pct_counts_ribo_low_outlier = adata.obs[adata.obs['ribo_outlier']==True]['pct_counts_ribo'].min()
            log1p_total_counts_high_outlier = adata.obs[adata.obs['outlier']==True]['log1p_total_counts'].max()
            log1p_total_counts_low_outlier = adata.obs[adata.obs['outlier']==True]['log1p_total_counts'].min()
            log1p_n_genes_by_counts_high_outlier = adata.obs[adata.obs['outlier']==True]['log1p_n_genes_by_counts'].max()
            log1p_n_genes_by_counts_low_outlier = adata.obs[adata.obs['outlier']==True]['log1p_n_genes_by_counts'].min()
            total_counts_high_outlier = adata.obs[adata.obs['outlier']==True]['total_counts'].max()
            total_counts_low_outlier = adata.obs[adata.obs['outlier']==True]['total_counts'].min()
            n_genes_by_counts_high_outlier = adata.obs[adata.obs['outlier']==True]['n_genes_by_counts'].max()
            n_genes_by_counts_low_outlier = adata.obs[adata.obs['outlier']==True]['n_genes_by_counts'].min()

            #Possible out of tissue spots
            # List of thresholds
            thresholds = [50, 100, 500, 1000, 5000, 10000, 20000]
            labels = [f'Below {t}' for t in thresholds] + ['Above 20000']
            # Create the new 'total_counts_category' based on the thresholds
            adata.obs['total_counts_thresholds'] = pd.cut(
                adata.obs['total_counts'],
                bins=[-1] + thresholds + [adata.obs['total_counts'].max()],
                labels=labels,
                right=True
            )

            # Save outlier information to a text file, including sample name
            with open(os.path.join(output_dir, 'outlier_value_counts.txt'), 'w') as f:
                f.write(f"Summary statistics for sample: {sample_name}\n")
                f.write(summary_stats.to_string() + "\n\n")
#                f.write(specific_stats.to_string() + "\n\n")
                f.write(f"Outlier value counts for sample: {sample_name}\n")
                f.write(str(adata.obs.outlier.value_counts()) + "\n\n")
                f.write(f"Mitochondrial outlier value counts for sample: {sample_name}\n")
                f.write(str(adata.obs.mt_outlier.value_counts()) + "\n\n")
                f.write(f"Ribosomal outlier value counts for sample: {sample_name}\n")
                f.write(str(adata.obs.ribo_outlier.value_counts()) + "\n\n")

                # Print high and low values for non-outliers for each metric
                f.write("High and Low Values for Non-Outliers:\n")
                f.write(f"pct_counts_mt - High: {pct_counts_mt_high}, Low: {pct_counts_mt_low}\n")
                f.write(f"pct_counts_ribo - High: {pct_counts_ribo_high}, Low: {pct_counts_ribo_low}\n")
                f.write(f"log1p_total_counts - High: {log1p_total_counts_high}, Low: {log1p_total_counts_low}\n")
                f.write(f"log1p_n_genes_by_counts - High: {log1p_n_genes_by_counts_high}, Low: {log1p_n_genes_by_counts_low}\n")
                f.write(f"total_counts - High: {total_counts_high}, Low: {total_counts_low}\n")
                f.write(f"n_genes_by_counts - High: {n_genes_by_counts_high}, Low: {n_genes_by_counts_low}\n")

                f.write("High and Low Values for Outliers:\n")
                f.write(f"pct_counts_mt - High: {pct_counts_mt_high_outlier}, Low: {pct_counts_mt_low_outlier}\n")
                f.write(f"pct_counts_ribo - High: {pct_counts_ribo_high_outlier}, Low: {pct_counts_ribo_low_outlier}\n")
                f.write(f"log1p_total_counts - High: {log1p_total_counts_high_outlier}, Low: {log1p_total_counts_low_outlier}\n")
                f.write(f"log1p_n_genes_by_counts - High: {log1p_n_genes_by_counts_high_outlier}, Low: {log1p_n_genes_by_counts_low_outlier}\n")
                f.write(f"total_counts - High: {total_counts_high_outlier}, Low: {total_counts_low_outlier}\n")
                f.write(f"n_genes_by_counts - High: {n_genes_by_counts_high_outlier}, Low: {n_genes_by_counts_low_outlier}\n")
                f.write(f"Spots for each threshold category:\n{adata.obs['total_counts_thresholds'].value_counts()}\n")


            #plot spatially the QC metrics
            plt.rcParams["figure.figsize"] = (8, 8)
            sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"], show=False)
            plt.savefig(os.path.join(output_dir, 'spatial_counts.png'), bbox_inches='tight')
            plt.close()
            sc.pl.spatial(adata, img_key="hires", color=["pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"], show=False)
            plt.savefig(os.path.join(output_dir, 'spatial_pct.png'), bbox_inches='tight')
            plt.close()
            
            sc.pl.spatial(adata, color='total_counts_thresholds', title='Total Counts Category')
            plt.savefig(os.path.join(output_dir, 'spatial_thresholds.png'), bbox_inches='tight')
            plt.close()


            #save adata
            adata.write(os.path.join(output_data_dir, "qc_metrics.h5ad"))
