import anndata as ad
import numpy as np
import h5py
import scanpy as sc
import pandas as pd

# ATAC DATA FORMATTING FOR PEAKS 

# DOWNLOAD ATAC gene by cell matrix into anndata structure
atac_path = 'data/atac/cerebrum-study/'

peak_data = ad.read_mtx(atac_path+'aggregate-section-8-peak-matrix.mtx')
peaknames = pd.read_csv(atac_path+'peaks-aggregate-section-8-peak-matrix.csv')['x']
metadata = pd.read_csv(atac_path+'section-8-metadata.csv')

# Extract chrom data from names
chroms = [s.split(':')[0] for s in peaknames]
chroms_start = [s.split(':')[1].split('-')[0] for s in peaknames]
chroms_end = [s.split(':')[1].split('-')[1] for s in peaknames]

peak_data.obs = pd.DataFrame({'Barcode':metadata['barcode']})
peak_data.var = pd.DataFrame({'chrom':chroms, 'chromStart':chroms_start,'chromEnd':chroms_end})


# Read metadata for all the samples 
files = list(("CEMBA180426_8B","CEMBA190711_8J","CEMBA180430_8B",
             "CEMBA190716_8E","CEMBA190711_8E","CEMBA190716_8J"))

all_metadata = pd.read_csv('data/atac/cerebrum-study/CEMBA-metadata.tsv',sep='\t')
metadata_section_8_from_all = all_metadata[all_metadata['Sample'].isin(files)].reset_index(drop=True)

# Download cell metadata and filter out non-included samples. 

# Filter out unannotated barcodes
filtered_peak_data = peak_data[peak_data.obs['Barcode'].isin(metadata_section_8_from_all['Barcode']),:]
filtered_peak_data.obs = filtered_peak_data.obs.reset_index(drop=True)

filtered_peak_data = filtered_peak_data[np.argsort(filtered_peak_data.obs['Barcode']),:]
filtered_peak_data.obs = filtered_peak_data.obs.reset_index(drop=True)

# Sort all metadata and append to adata
filtered_peak_data = filtered_peak_data[np.argsort(filtered_peak_data.obs['Barcode']),:]
metadata_section_8_from_all = metadata_section_8_from_all.sort_values(by=['Barcode']).reset_index(drop=True)

if all(filtered_peak_data.obs['Barcode'] == metadata_section_8_from_all['Barcode']):
    filtered_peak_data.obs = metadata_section_8_from_all
    print('Properly sorted and integrated')
    
# ATAC-seq data in anndata format ready for analyzing

filtered_peak_data.write_h5ad(atac_path+'aggregate-section-8-valid-peak-count-data.h5ad')
