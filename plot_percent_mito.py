import scanpy as sc
import scvelo as scv
import numpy as np

adata = scv.read('data/full_joint.h5ad', cache=True)
scv.utils.show_proportions(adata)
scv.pl.proportions(adata, groupby='type')

percent_mito_lowerbd = [0.5, 1, 2, 5]
for i in percent_mito_lowerbd:
	adata.obs['percent.mito.>{}'.format(i)] = (adata.obs['percent.mito'] > i).astype(int)
	print(round(np.where(adata.obs['percent.mito.>{}'.format(i)] == 0)[0].size / adata.shape[0], 5))
	sc.pl.embedding(adata, 'umap', color=['percent.mito.>{}'.format(i)], save='_full_pct_mito_{}.png'.format(i))
