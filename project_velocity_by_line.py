import scvelo as scv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('subset', help='the subset of cells where we are performing the analysis, i.e. E1, E1col1, or full')
parser.add_argument('velocity_mode', nargs='?', default='dynamical', help='deterministic, stochastic or dynamical')
parser.add_argument('plot_mode', nargs='?', default='stream', help='stream or arrow')
parser.add_argument('color_mode', nargs='?', default='type', help='type or day')
parser.add_argument('tkey', nargs='?', type=int, default=0, help='whether to use days as `tkey` when computing velocity graph')
args = parser.parse_args()


lines = [
    '18489', '18499', '18505', '18508', '18511', '18517', '18520',
    '18855', '18858', '18870', '18907', '18912', '19093', '19108', 
    '19127', '19159', '19190', '19193', '19209'
]

data_path = 'data/{}_{}_velocity.h5ad'.format(args.subset, args.velocity_mode)
adata = scv.read(data_path, cache=True)
scv.tl.velocity_embedding(adata, basis='umap')

for i, line in enumerate(lines):
    subset = adata[adata.obs['individual'] == line]
    if args.plot_mode == 'stream':
        scv.pl.velocity_embedding_stream(subset, basis='umap', dpi=100, title=line, color=args.color_mode,
            save='{}_{}_{}_{}_{}_tkey{}.png'.format(args.subset, args.plot_mode, args.color_mode, args.velocity_mode, line, args.tkey))
    elif args.plot_mode == 'arrow':
        scv.pl.velocity_embedding(subset, basis='umap', dpi=100, title=line, color=args.color_mode,
            save='{}_{}_{}_{}_{}_tkey{}.png'.format(args.subset, args.plot_mode, args.color_mode, args.velocity_mode, line, args.tkey))
