import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc
# RNA DATA FORMATTING
rna_data_path = 'data/rna/cell2location/'
samples = ['5705STDY8058280','5705STDY8058281','5705STDY8058282','5705STDY8058283','5705STDY8058284','5705STDY8058285']

adata_list = []
for sample in samples:
    adata = sc.read_10x_h5(rna_data_path + sample+'_filtered_feature_bc_matrix.h5')
    adata.obs['Cell ID'] = [sample +'_'+ idx for idx in adata.obs.index]
    adata.var_names_make_unique()
    adata_list.append(adata)

adata_rna = ad.concat(adata_list)
metadata = pd.read_csv(rna_data_path+'cell_annotation.csv')



# Merge & Filter 
adata_rna_filtered = adata_rna[adata_rna.obs['Cell ID'].isin(metadata_ordered['Cell ID']),:]
metadata_ordered = pd.merge(adata_rna.obs,metadata,'inner',on='Cell ID',)


if all(adata_rna_filtered.obs['Cell ID'].reset_index(drop=True) == metadata_ordered['Cell ID']):
    adata_rna_filtered.obs = metadata_ordered
    

adata_rna_filtered.write(rna_data_path+'aggregate_cell2location_scrna_samples_annotated.h5ad')