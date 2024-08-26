#Code to run alignment of both spatial transcriptomics datasets using PASTE: to use activate conda environment "squidpy-env"


import math
import time
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import style
import paste as pst
import os
import squidpy
import anndata as ad

from skimage import io



data_path = '/home/quiquevb/Desktop/Data/Yadav'

for folder in os.listdir(data_path):
    if folder == 'Sample1':
        sample_path = os.path.join(data_path, folder)
        counts_file = os.path.join(sample_path, 'GSM6919905_V160-CGND-HRA-02636-C_filtered_feature_bc_matrix.h5')
    elif folder == 'Sample2_original':
        sample_path2 = os.path.join(data_path, folder)
        counts_file2 = os.path.join(sample_path2, 'GSM6919906_V160-CGND-HRA-02636-D_filtered_feature_bc_matrix.h5')

adata_1 = squidpy.read.visium(path=sample_path, counts_file=counts_file)
adata_2 = squidpy.read.visium(path=sample_path2, counts_file=counts_file2)

adata_1.var_names_make_unique()
adata_2.var_names_make_unique()

start = time.time()
pi12 = pst.pairwise_align(adata_1, adata_2)

pis = [pi12]
slices = [adata_1, adata_2]
new_slices, angles, translation  = pst.stack_slices_pairwise(slices, pis, output_params=True, matrix=True)


#The following plot shows alignment
slice_colors = ['#e41a1c','#377eb8']

tissue_image = "/home/quiquevb/Desktop/Data/Yadav/Sample1/GSM6919905_V160-CGND-HRA-02636-C_tissue_hires_image.png"
plt.figure(figsize=(7,7))
for i in range(len(new_slices)):
    pst.plot_slice(new_slices[i],slice_colors[i])
plt.legend(handles=[mpatches.Patch(color=slice_colors[0], label='1'),mpatches.Patch(color=slice_colors[1], label='2')])

plt.gca().invert_yaxis()
plt.axis('off')
plt.show()


# Assuming adata1 and adata2 are your two AnnData objects
adata1 = new_slices[0]  # The first AnnData object
adata2 = new_slices[1]  # The second AnnData object

# Add a new column to `obs` indicating the slice number
adata1.obs['slice'] = 'slice_1'
adata2.obs['slice'] = 'slice_2'

# Concatenate the two AnnData objects
combined_adata = ad.concat([adata1, adata2])

combined_adata.write_h5ad(os.path.join(data_path, "combined_slices.h5ad"))

plt.rcParams["figure.figsize"] = (3, 3)
combined_adata.obsm['spatial'][:, 1] = -1*combined_adata.obsm['spatial'][:, 1]
combined_adata.obs['slice'].replace({'S1':'slice_1', 'S3':'slice_2'}, inplace=True)
ax = sc.pl.embedding(combined_adata, basis='spatial',
                color='slice',
                show=False)
ax.set_title('Aligned image')
#ax.axis('off')
