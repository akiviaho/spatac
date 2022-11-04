from itertools import chain
import time

import anndata as ad
import itertools
import pandas as pd
import scanpy as sc
import seaborn as sns
import numpy as np
from scipy import sparse

date = 20220929
sp_data_folder ='synthetic-spatial/'
filter_cells = True

# Load data
synthetic_cell_counts = pd.read_csv('synthetic-spatial/cell_count_df_20220929.csv',index_col=0)
synthetic_cell_counts = synthetic_cell_counts.astype(np.int16)
atac = ad.read_h5ad('data/share-seq/mouse-brain/share-seq-mouse-brain-atac-data.h5ad')

# CHANGE CELL TYPE ANNOTATIONS INTO A BROADER SPECTRUM
if filter_cells:
    cells_to_keep = ['EN','IN','A1.E1','OG1']
    cell_idxs = np.where(atac.obs['celltype'].str.contains('|'.join(cells_to_keep)))[0]
    atac = atac[cell_idxs]

    aggr_celltypes = atac.obs['celltype'].astype(str)
    aggr_celltypes[aggr_celltypes.str.contains('EN')] = 'EN'
    aggr_celltypes[aggr_celltypes.str.contains('IN')] = 'IN'
    atac.obs['broad_celltype'] = aggr_celltypes


# Create empty numpy matrix for peak data
spots_by_peaks = np.empty((len(synthetic_cell_counts),len(atac.var)),dtype=np.int16)
celltypes = synthetic_cell_counts.columns

for row_idx, row in synthetic_cell_counts.reset_index(drop=True).iterrows():
    select_indices = np.empty((0))

    for col_idx,cell_count in enumerate(row): # Iterates over columns in each row
        celltype = celltypes[col_idx]
        celltype_indices = np.where(atac.obs.broad_celltype==celltype)[0] # Finds ATAC data with cell types

        # Random select from matching indices according to number of cells present
        select_indices = np.append(select_indices,np.random.choice(celltype_indices,cell_count)) 

    # Subsets ATAC by selected index and produces a sum of those cells' peaks, appends to 
    # correct synthetic spot
    spots_by_peaks[row_idx,:] = np.sum(atac.X[select_indices,:].toarray(),axis=0,dtype=np.int16)

    
# Save peak counts to 
synthetic_spatial_atac_data = ad.AnnData(X=sparse.csr_matrix(spots_by_peaks),
                                  obs=synthetic_cell_counts,
                                  var=atac.var)

synthetic_spatial_atac_data.write(f'{sp_data_folder}synth_spatial_atac_data_from_cell_counts_{date}_sparse.h5ad')