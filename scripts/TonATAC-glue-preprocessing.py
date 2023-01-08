import anndata as ad
import networkx as nx
import scanpy as sc
import pandas as pd
import numpy as np
import scglue
from matplotlib import rcParams

import logging
import os
os.chdir('/lustre/scratch/kiviaho/spatac/integrations/tonsilatlas/actual_scatac_to_spatial')
date = '20221215'
multiple_experiments = False

rna_name = 'actual_spatial_rna'
atac_name = 'simulated_spots_atac'

rna = ad.read_h5ad(rna_name+'.h5ad')
atac = ad.read_h5ad(atac_name + '.h5ad')

# Copy raw counts into X, only if 
# For single cell
#rna.X = rna.raw.X.copy()
#atac.X = atac.raw.X.copy()

# For simulated data
rna.layers['counts'] = rna.X.copy()
atac.layers['counts'] = atac.X.copy()

scglue.data.get_gene_annotation(
    rna, gtf='../../../gencode.v42.annotation.gtf.gz',
    gtf_by="gene_name"
)
logging.warning(' Annotation built')

# Drop unannotated genes:
rna = rna[:,~rna.var.index.duplicated(keep='first')]
rna = rna[:,rna.var.dropna(subset=['chrom','chromStart','chromEnd']).index]
sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3",span=1)
sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.scale(rna)
sc.tl.pca(rna, n_comps=100, svd_solver="auto")

logging.warning(' RNA preprocessing done')

# Drop peak locations with zero peaks
atac = atac[:,~(atac.X.sum(axis=0)==0)]
atac = atac[~(atac.X.sum(axis=1)==0),:] ## drop sample locations with zero counts

# Rename genomic coordinates to match GLUE scheme
atac.var.rename(columns={'seqnames':'chrom','start':'chromStart','end':'chromEnd'},inplace=True) # just rename existing columns !

# Subset the anndatas for graph building, since it doesn't need data points but still uses a shit ton of memory
sub_rna = rna[:5]
sub_atac = atac[:5]
sub_atac.var.rename(columns={'seqnames':'chrom','start':'chromStart','end':'chromEnd'},inplace=True) # just rename existing columns !

# Moved guidance graph construction to BEFORE lsi calculation, since it should be using HVGs (which it gets from the RNA)
logging.warning(' Constructing guidance graph...')
guidance = scglue.genomics.rna_anchored_guidance_graph(sub_rna, sub_atac)
scglue.graph.check_graph(guidance, [sub_rna, sub_rna])

atac.var = sub_atac.var.copy()

logging.warning(' Checking if ATAC contains highly variable information:')
atac.var.head()

logging.warning(' Calculating ATAC lsi...')
scglue.data.lsi(atac,n_components=100) # This takes a long time, but not so long when hvf's have been determined (from RNA)
logging.warning(' ATAC lsi done')


logging.warning(' Writing ready ATAC file...')
atac.write('preprocessed_'+atac_name+'_'+date+".h5ad", compression="gzip")

logging.warning(' Writing graph...')
nx.write_graphml(guidance, "guidance_graph_"+atac_name+"_"+rna_name+"_"+date+".graphml.gz")

logging.warning(' Writing ready RNA file...')
try:
    rna.write("preprocessed_"+rna_name+"_"+date+".h5ad", compression="gzip")
except:
    rna.var = rna.var.drop(columns='artif_dupl')
    rna.write("preprocessed_"+rna_name+"_"+date+".h5ad", compression="gzip")
    