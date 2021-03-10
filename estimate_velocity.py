'''
Preprocess the spliced/unspliced data and estimate velocity
'''

import scanpy as sc
import scvelo as scv
import scachepy
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('subset', help='the subset of cells where we are performing the analysis, i.e. E1, E1col1, or full')
parser.add_argument('mode', nargs='?', default='dynamical', help='deterministic, stochastic or dynamical')
args = parser.parse_args()

if args.mode not in ['deterministic', 'stochastic', 'dynamical']:
	raise ValueError('{} is not a velocity mode'.format(args.mode))

input_path = 'data/{}_joint.h5ad'.format(args.subset)
output_path = 'data/{}_{}_velocity.h5ad'.format(args.subset, args.mode)
cache_dir = 'cached_files/{}'.format(args.subset)

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo')

adata = scv.read(input_path, cache=True)
scv.utils.show_proportions(adata)

if args.mode == 'dynamical':
	c = scachepy.Cache(cache_dir)
	c.tl.recover_dynamics(adata, force=False, n_jobs=16)
	scv.tl.velocity(adata, mode='dynamical')
else:
	scv.tl.velocity(adata, mode=args.mode)

scv.tl.velocity_graph(adata)

adata.write(output_path)



