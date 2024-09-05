import os
import torch
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp
import matplotlib.pyplot as plt
from GraphST import GraphST

# Run device, by default, the package is implemented on 'cpu'. We recommend using GPU.
device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')

# the location of R, which is necessary for mclust algorithm. Please replace it with local R installation path
os.environ['R_HOME'] = '/home/quiquevb/.conda/envs/graphst/lib/R'

#read data
file_fold = '/storage/gge/Quique/Yadav/spatial_data'
output_dir = '/storage/gge/Quique/Yadav/spatial_data/Output'

adata = sc.read_h5ad(file_fold + '/combined_slices.h5ad')
adata.var_names_make_unique()

plt.rcParams["figure.figsize"] = (3, 3)
adata.obsm['spatial'][:, 1] = -1*adata.obsm['spatial'][:, 1]
adata.obs['slice'].replace({'S1':'slice_1', 'S3':'slice_2'}, inplace=True)
ax = sc.pl.embedding(adata, basis='spatial',
                color='slice',
                show=False)
ax.set_title('Aligned image')

### Plotting UMAP before batch effect correction
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)

sc.pp.neighbors(adata, use_rep='X_pca', n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.pl.umap(adata, color='slice', title='Uncorrected',
                  show=False)

plt.savefig(os.path.join(output_dir, "UMAP_no_batch_correction.png"))

#ax.axis('off')


# define model
model = GraphST.GraphST(adata, device=device)

# run model
adata = model.train()


# set radius to specify the number of neighbors considered during refinement
radius = 50
n_clusters = 7

tool = 'mclust' # mclust, leiden, and louvain

# clustering
from GraphST.utils import clustering

if tool == 'mclust':
   clustering(adata, n_clusters, radius=radius, method=tool, refinement=True) # For DLPFC dataset, we use optional refinement step.
elif tool in ['leiden', 'louvain']:
   clustering(adata, n_clusters, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=False)

adata.write(os.path.join(file_fold, "combined_slices_clustering.h5ad"))
