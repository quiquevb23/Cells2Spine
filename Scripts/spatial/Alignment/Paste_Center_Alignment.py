#Script to center align all 4 slices using PASTE

import os
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import paste as pst

# Define directories
base_dir = '/storage/gge/Quique/Cells2SpineData/Pilot/spatial/'
output_base_dir = '/home/quiquevb/Cells2Spine/Cells2Spine/Outputs/spatial'
slices = []

# Create output directories
output_dir = os.path.join(output_base_dir, 'Align')
os.makedirs(output_dir, exist_ok=True)

# Loop over all subdirectories
for file_name in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, file_name)
    for file in os.listdir(sample_path):
        if file.endswith('h5ad'):
            sample_name = file.replace('h5ad', '')
            counts_file = os.path.join(sample_path, file)

            # Load the data
            f'adata_{sample_name}' = sc.read(counts_file)
            slices.append(f'adata{sample_name}')

lmbda = len(slices)*[1/len(slices)]

pst.filter_for_common_genes(slices)

b = []
for i in range(len(slices)):
    b.append(pst.match_spots_using_spatial_heuristic(slices[0].X, slices[i].X))

start = time.time()

center_slice, pis = pst.center_align(initial_slice, slices, lmbda, random_seed = 5, pis_init = b)

print('Runtime: ' + str(time.time() - start))

#now plot results

center, new_slices = pst.stack_slices_center(center_slice, slices, pis)

center_color = 'orange'
slices_colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3']

plt.figure(figsize=(7,7))
pst.plot_slice(center,center_color,s=400)
for i in range(len(new_slices)):
    pst.plot_slice(new_slices[i],slices_colors[i],s=400)

plt.legend(handles=[mpatches.Patch(color=slices_colors[0], label='1'),mpatches.Patch(color=slices_colors[1], label='2'),mpatches.Patch(color=slices_colors[2], label='3'),mpatches.Patch(color=slices_colors[3], label='4')])
plt.gca().invert_yaxis()
plt.axis('off')
plt.show()
