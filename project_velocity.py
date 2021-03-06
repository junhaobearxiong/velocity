import scanpy as sc
import scvelo as scv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('subset', help='the subset of cells where we are performing the analysis, i.e. E1, E1col1, or full')
parser.add_argument('velocity_mode', nargs='?', default='dynamical', help='deterministic, stochastic or dynamical')
parser.add_argument('plot_mode', nargs='?', default='stream')
parser.add_argument('color_mode', nargs='?', default='type')
parser.add_argument('tkey', nargs='?', type=int, default=0, help='whether to use days as `tkey` when computing velocity graph')
args = parser.parse_args()

if not args.tkey:
    data_path = 'data/{}_{}_velocity.h5ad'.format(args.subset, args.velocity_mode)
else:
    data_path = 'data/{}_{}_velocity_tkey.h5ad'.format(args.subset, args.velocity_mode)

adata = sc.read(data_path)

if args.plot_mode == 'stream':
    scv.pl.velocity_embedding_stream(adata, basis='umap', legend_fontsize=12, dpi=200, title=args.subset, color=args.color_mode, 
        save='{}_{}_{}_{}_tkey{}.png'.format(args.subset, args.plot_mode, args.color_mode, args.velocity_mode, args.tkey))

elif args.plot_mode == 'arrow':
    scv.pl.velocity_embedding(adata, basis='umap', dpi=200, density=0.75, title=args.subset, color=args.color_mode, 
        save='{}_{}_{}_{}_tkey{}.png'.format(args.subset, args.plot_mode, args.color_mode, args.velocity_mode, args.tkey))
