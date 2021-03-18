import scanpy as sc
import scvelo as scv
import cellrank as cr
import scachepy
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('subset', help='the subset of cells where we are performing the analysis, i.e. E1, E1col1, or full')
parser.add_argument('mode', nargs='?', default='dynamical', help='deterministic, stochastic or dynamical')
args = parser.parse_args()

if args.mode not in ['deterministic', 'stochastic', 'dynamical']:
	raise ValueError('{} is not a velocity mode'.format(args.mode))

input_path = 'data/{}_{}_velocity.h5ad'.format(args.subset, args.mode)
output_path = 'data/{}_{}_cellrank.h5ad'.format(args.subset, args.mode)
cache_dir = 'cached_files/{}_{}'.format(args.subset, args.mode)

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo')
cr.settings.verbosity = 2

adata = scv.read(input_path, cache=True)
scv.utils.show_proportions(adata)

# plot terminal and initial states 
cr.tl.terminal_states(adata, cluster_key='type', weight_connectivities=0.2)
cr.tl.initial_states(adata, cluster_key='type')
cr.pl.terminal_states(adata, dpi=200, save='{}_{}_terminal_states.png'.format(args.subset, args.mode))
cr.pl.initial_states(adata, dpi=200, save='{}_{}_initial_states.png'.format(args.subset, args.mode))

# plot lineages, i.e. fate probabilities of each cell
cr.tl.lineages(adata)
cr.pl.lineages(adata, dpi=200, same_plot=False, save='{}_{}_lineages_diffplot.png'.format(args.subset, args.mode))
cr.pl.lineages(adata, dpi=200, same_plot=True, save='{}_{}_lineages_sameplot.png'.format(args.subset, args.mode))

adata.write(output_path)
