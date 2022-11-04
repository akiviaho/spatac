from itertools import chain

import anndata as ad
import torch
import itertools
import networkx as nx
import pandas as pd
import scanpy as sc
import scglue
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

multiple_experiments = False
exp = '20220929'

if __name__ == "__main__":    
    
    if torch.cuda.is_available():
        print('GPU support available. Proceeding...')
    else:
        print('GPU support not available. Exiting...')
    
    print('Downloading atac data & guidance graph...')
    
    atac = ad.read_h5ad('preprocessed_shareseq_atac_data_20220929.h5ad')
    guidance = nx.read_graphml('guidance_synth_spatial_atac_shareseq_20220929.graphml.gz')
    rna = ad.read_h5ad('preprocessed_synthetic_spatial_data_20220929.h5ad')
    
#     if multiple_experiments:
#         experiment_ids = ['exper'+ str(n) for n in list(range(0,10))]
#         #experiment_ids = ['all_experiments']
#         for exp in experiment_ids:

#     rna = ad.read_h5ad('data/preprocessed_synthetic_spatial_'+exp+'.h5ad')

    ####### MODEL CONFIGURATION AND TRAINING ########
    print('Configuring the model...')

    hostdevice = '/device:gpu:0'

    scglue.models.configure_dataset(
        rna, "NB", use_highly_variable=True,
        use_layer="counts", use_rep="X_pca"
    )

    scglue.models.configure_dataset(
        atac, "NB", use_highly_variable=True,
        use_rep="X_lsi"
    )

    guidance_hvf = guidance.subgraph(chain(
        rna.var.query("highly_variable").index,
        atac.var.query("highly_variable").index
    )).copy()

    print('Training the model...')
    glue = scglue.models.fit_SCGLUE(
        {"rna": rna, "atac": atac}, guidance_hvf,
        fit_kws={"directory": "glue"}
    )

    ## SAVE GLUE MODEL ##

    print('Saving the model...')
    glue.save("results/glue-model-shareseq-synthetic-"+exp+'.dill')
    print('Successfully saved!')

    dx = scglue.models.integration_consistency(
        glue, {"rna": rna, "atac": atac}, guidance_hvf
    )

    rna.obsm["X_glue"] = glue.encode_data("rna", rna)
    atac.obsm["X_glue"] = glue.encode_data("atac", atac)


    # Lineplot of consistency score
    lineplot = sns.lineplot(x="n_meta", y="consistency", data=dx).axhline(y=0.05, c="darkred", ls="--")
    fig = lineplot.get_figure()
    fig.savefig('results/glue-model-consistency-score-'+exp+'.png',dpi=200)

    feature_embeddings = glue.encode_graph(guidance_hvf)
    feature_embeddings = pd.DataFrame(feature_embeddings, index=glue.vertices)

    rna.varm["X_glue"] = feature_embeddings.reindex(rna.var_names).to_numpy()
    atac.varm["X_glue"] = feature_embeddings.reindex(atac.var_names).to_numpy()


    # SAVE THE EMBEDDINGS
    print('Saving data with embeddings')
    rna.write('results/rna-synthetic-spatial-shareseq-'+exp+'-with-glue-embeddings.h5ad', compression="gzip")
    atac.write('results/atac-synthetic-spatial-shareseq-'+exp+'-with-glue-embeddings.h5ad', compression="gzip")
    nx.write_graphml(guidance_hvf, "results/synthetic-spatial-shareseq-guidance-hvf-"+exp+".graphml.gz")

    print('Plotting UMAPs')

    combined = ad.concat([rna, atac])
    combined.obs['domain'] = np.concatenate([np.repeat('scRNA-seq',len(rna)),np.repeat('scATAC-seq',len(atac))])

    sc.pp.neighbors(combined, use_rep="X_glue", metric="cosine")
    sc.tl.umap(combined)
    sc.set_figure_params(scanpy=True, dpi=100, dpi_save=150, vector_friendly=True, fontsize=14, figsize=(12,12))
    sc.pl.umap(combined, color=["domain"], wspace=0.65,save='scglue-shared-embedding-'+exp+'.png')


