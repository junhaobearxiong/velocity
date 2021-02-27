import scanpy as sc
import scvelo as scv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('subset', help='the subset of cells where we are performing the analysis, i.e. E1, E1col1, or full')
parser.add_argument('tkey', type=int, nargs='?', default=0, help='whether to use days as `tkey` when computing velocity graph')
args = parser.parse_args()

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo')

if not args.tkey:
    data_path = 'data/{}_velocity.h5ad'.format(args.subset)
else:
    data_path = 'data/{}_velocity_tkey.h5ad'.format(args.subset)

adata = sc.read(data_path)
scv.tl.terminal_states(adata)
adata.write(data_path)