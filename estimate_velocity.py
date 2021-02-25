'''
Preprocess the spliced/unspliced data and estimate velocity
'''

import scanpy as sc
import scvelo as scv
import scachepy
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('subset', help='the subset of cells where we are performing the analysis, i.e. E1, E1col1, or full')
args = parser.parse_args()

input_path = 'data/{}_joint.h5ad'.format(args.subset)
output_path = 'data/{}_velocity.h5ad'.format(args.subset)
cache_dir = 'cached_files/{}'.format(args.subset)

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo')

adata = scv.read(input_path, cache=True)
scv.utils.show_proportions(adata)

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

c = scachepy.Cache(cache_dir)
c.tl.recover_dynamics(adata, force=False, n_jobs=16)

scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
sc.tl.umap(adata)

adata.write(output_path)



