import anndata as ad
import numpy as np
import h5py
import scanpy as sc
import pandas as pd

# DOWNLOAD ATAC gene by cell matrix into anndata structure
atac_path = 'data/atac/cerebrum-study/'

gene = ad.read_mtx(atac_path+'aggregate-section-8-gene-matrix.mtx')
genenames = pd.read_csv(atac_path+'genes-aggregate-section-8-gene-matrix.csv')['x']
metadata = pd.read_csv(atac_path+'section-8-metadata.csv')

gene.obs = pd.DataFrame({'Barcode':metadata['barcode']})
gene.var = pd.DataFrame({'Gene':genenames})

# Download cell metadata and filter out non-included samples. 
files = list(("CEMBA180426_8B","CEMBA190711_8J","CEMBA180430_8B",
             "CEMBA190716_8E","CEMBA190711_8E","CEMBA190716_8J"))

all_metadata = pd.read_csv('data/atac/cerebrum-study/CEMBA-metadata.tsv',sep='\t')
metadata_section_8_from_all = all_metadata[all_metadata['Sample'].isin(files)].reset_index(drop=True)

# Filter out unannotated barcodes
filtered_gene = gene[gene.obs['Barcode'].isin(metadata_section_8_from_all['Barcode']),:]
filtered_gene.obs = filtered_gene.obs.reset_index(drop=True)

filtered_gene = filtered_gene[np.argsort(filtered_gene.obs['Barcode']),:]
filtered_gene.obs = filtered_gene.obs.reset_index(drop=True)

# Sort all metadata and append to adata
filtered_gene = filtered_gene[np.argsort(filtered_gene.obs['Barcode']),:]
metadata_section_8_from_all = metadata_section_8_from_all.sort_values(by=['Barcode']).reset_index(drop=True)

if all(filtered_gene.obs['Barcode'] == metadata_section_8_from_all['Barcode']):
    filtered_gene.obs = metadata_section_8_from_all
    print('Properly sorted and integrated')

# ATAC-seq data in anndata format ready for analyzing
filtered_gene.var.index = filtered_gene.var['Gene']
# Write to correct path
filtered_gene.write_h5ad(atac_path+'aggregate-section-8-valid-cells.h5ad')