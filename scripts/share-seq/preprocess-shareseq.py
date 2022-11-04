import anndata as ad
import pandas as pd
import numpy as np

# DOWNLOAD RNA DATA FROM SHARE-SEQ EXPERIMENT
rna = pd.read_csv('data/share-seq/mouse-brain/GSM4156610_brain.rna.counts.txt',sep='\t')
rna = rna.set_index('gene')

# Build into anndata
rna = ad.AnnData(X=np.array(rna).T,
                       var=pd.DataFrame({'Gene':np.asarray(rna.index)}),
                      obs=pd.DataFrame({'Barcode':np.asarray([bc.replace(',','.') for bc in rna.columns])}))

# DOWNLOAD ATAC DATA FROM SHARE-SEQ EXPERIMENT
peaks = pd.read_csv('data/share-seq/mouse-brain/GSM4156599_brain.peaks.bed',sep='\t',header=None)
peaks.rename(columns={0:'chrom',1:'chromBegin',2:'chromEnd'},inplace=True)
celltypes = pd.read_csv('data/share-seq/mouse-brain/GSM4156599_brain_celltype.txt',sep='\t')
counts = ad.read_mtx('data/share-seq/mouse-brain/GSM4156599_brain.counts.txt').T

# Build the ATAC anndata structure
atac = ad.AnnData(X=counts)
atac.obs = celltypes
atac.var = peaks

rna.write('data/share-seq/mouse-brain/share-seq-mouse-brain-rna-data.h5ad', compression="gzip")
atac.write('data/share-seq/mouse-brain/share-seq-mouse-brain-atac-data.h5ad', compression="gzip")