import scanpy as sc
import scvelo as scv

adata = scv.read('data/E1col1_joint.h5ad', cache=True)
scv.utils.show_proportions(adata)
scv.pl.proportions(adata, groupby='type')

percent_mito_lowerbd = [0.1, 0.5, 1, 5]
for i in percent_mito_lowerbd:
	adata.obs['percent.mito.>{}'.format(i)] = (adata.obs['percent.mito'] > i).astype(int)
	sc.pl.embedding(adata, 'umap', color=['percent.mito.>1'], dpi=150, save='_E1col1_pctmito_{}.png'.format(i))
