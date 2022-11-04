import anndata as ad
import networkx as nx
import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.decomposition import TruncatedSVD
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize
import scglue

n_dims = 100
n_barcodes = 'all'

def filter_rna_by_class(adata):
# Insert broad class labels into rna data
    rna_annot = adata.obs['annotation_1']
    rna_broad_class = np.repeat('unannotated',len(rna_annot))

    for i,a in enumerate(rna_annot):
        if 'Astro' in a:
            rna_broad_class[i] = 'Astrocyte'
        if 'Ext' in a:
            rna_broad_class[i] = 'Excitatory'
#        if 'Inh' in a:
#            rna_broad_class[i] = 'Inhibitory'
        if 'Oligo' in a:
            rna_broad_class[i] = 'Oligodendrocyte'
        if 'OPC' in a:
            rna_broad_class[i] = 'OPC'
        if 'Micro' in a:
            rna_broad_class[i] = 'Microglia'

    adata.obs['Broad Class'] = rna_broad_class

    adata = adata[(adata.obs['Broad Class']) != 'unannotated',:]
    return adata

def filter_atac_by_class(adata):
    # Insert broad class labels into atac data
    atac_broad_class = pd.read_csv('data/atac/cerebrum-study/cell-type-annotations.tsv',sep='\t',header=0)
    adata_major_type = adata.obs['MajorType']
    atac_meta = pd.merge(adata_major_type,atac_broad_class,how='left')
    adata.obs['Broad Class'] = np.asarray(atac_meta['Broad Class'])
    adata = adata[(adata.obs['Broad Class'] != 'unannotated'),:]
    return adata

def normalize_and_pca(adata,n_comp=100):
    
    print('Normalizing rna data')
    adata.layers["counts"] = adata.X.copy()
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3")
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.scale(adata)
    
    # Subset to HVG
    hvg_arr = np.asarray(adata[:,adata.var['highly_variable'] == True].X.T)
    hvg_arr = np.nan_to_num(hvg_arr)
    pca = PCA(n_components=n_comp)
    
    print('Calculating PCs for rna data')
    pca_arr = pca.fit(hvg_arr)
    
    adata.obsm['X_pca'] = pca_arr.components_.T
        
    print('PCA done!')
    return adata


def idf_normalize_and_lsi(adata,n_comp=100):
    #sc.pp.highly_variable_genes(adata, n_top_genes=10000, flavor="seurat_v3")
    
    print('TFIDF normalizing...')
    arr = adata.X.copy()
    print('Calculating idf')
    idf = arr.shape[0] / arr.sum(axis=0)
    print('Calculating tf')
    tf = arr.multiply(1 / arr.sum(axis=1))
    arr = tf.multiply(idf)
    print('Normalizing...')
    arr = normalize(arr, norm="l2")
    arr = np.log1p(arr * 1e4)
    adata.layers['normalized'] = arr
    
    # Subset to HVG and transform to dense
    print('Calculating variable genes')
    #hvg_arr = np.asarray(adata[:,adata.var['highly_variable'] == True].layers['normalized'].T)
    hvg_arr = np.nan_to_num(np.asarray(adata.layers['normalized'].T))
    
    print('Doing Singular Value Decomposition')
    svd = TruncatedSVD(n_components=n_comp,n_iter=30)
    svd_arr = svd.fit(hvg_arr)
    
    adata.obsm['X_lsi'] = svd_arr.components_.T
    
    return adata
    
if __name__ == "__main__":    

    rna = ad.read_h5ad('data/rna/cell2location/aggregate_cell2location_scrna_samples_annotated.h5ad')

    # Annotation downloaded from 
    # https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/gencode.vM30.annotation.gtf.gz
    scglue.data.get_gene_annotation(
        rna, gtf="gencode.vM30.annotation.gtf.gz",
        gtf_by="gene_name"
    )
    
    # Remove genes that did not have chromosome information or are mitochondrial
    rna = rna[:,np.where(rna.var['chrom'].notna())[0]]
    rna = rna[:,rna.var['chrom'] != 'chrM']
    
    rna = filter_rna_by_class(rna)
    rna = normalize_and_pca(rna)

    atac = ad.read_h5ad('data/atac/cerebrum-study/aggregate-section-8-valid-peak-count-data.h5ad')
    atac = filter_atac_by_class(atac)
    atac = idf_normalize_and_lsi(atac)

    #### INTEGRATE ATAC & RNA #####
    print('Calculating guidance graph')
    guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
    scglue.graph.check_graph(guidance, [rna, atac])

    rna.write(str(n_dims)+'_dimensions'+"_preprocessed-c2l-rna-data.h5ad", compression="gzip")
    atac.write(str(n_dims)+'_dimensions'+"_preprocessed-cereb-atac-data.h5ad", compression="gzip")
    nx.write_graphml(guidance, "guidance-rna-c2l-atac-cereb.graphml.gz")

