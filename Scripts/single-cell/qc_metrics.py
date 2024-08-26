import scanpy as sc
import anndata as ad
import pandas as pd


base_dir = '/storage/gge/Quique/Cells2SpineData/Pilot/single-cell/matrices/'
output_base_dir = '/home/quiquevb/Cells2Spine/Cells2Spine/Outputs/single-cell'

# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name)
    for file in os.listdir(sample_path):
        if file.endswith('filtered_feature_bc_matrix.h5'):
            output_dir = os.path.join(output_base_dir, sample_name, 'QC_metrics')
            sample_name = file.replace('_filtered_feature_bc_matrix.h5', '')
            adata = sc.read(os.path.join(sample_path, file)
            adata.layers['counts'] = adata.X #save raw counts
            adata.var['mt'] = adata.var_names.str.startswith('mt-')
            sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'],percent_top=None, log1p=False, inplace=True)
            #QC metrics are calculated with unnormalized data and stored, so if we want to use them later they are unnormalized

            sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],
              jitter=0.4, multi_panel=True, show=False)
            plt.savefig(os.path.join(output_dir, 'violin.png'))


            adata = adata[adata.obs['pct_counts_mt'] < 10, :].copy()
            adata = adata[adata.obs['n_genes_by_counts'] >500, :].copy()
            adata = adata[adata.obs['total_counts'] < 5000, :].copy()
            sc.pp.filter_genes(adata, min_cells=1) #filter out genes not expressed in any cell
            #sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'],percent_top=None, log1p=False, inplace=True)
            #we need to eliminate mito genes
            #now we need to normalize either with this or with seuratrecipe, select highly variable genes and then scale the data

            #we will not remove mitochondrial genes to better cluster
            sc.pp.normalize_total(adata, target_sum=10000)
            sc.pp.log1p(adata)
            #adata = adata[:, ~adata.var['mt']
            #update adata x from numpy array to scipy csr sparse matrix to avoid bug
            #adata.X = scipy.sparse.csr_matrix(adata.X)

            sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=5000)

            adata.raw = adata
            adata = adata[:, adata.var.highly_variable]

            #We will not regress out the effects on pct_counts_mt because the droputs get linearized and that is abnormal

            sc.pp.regress_out(adata, ["pct_counts_mt"])
            sc.pp.scale(adata)

            sc.tl.pca(adata, svd_solver="arpack")
            sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True,save='/raw_and_normalized_variance_pca.png')
            
