import scanpy as sc
import scvelo as scv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('subset', help='the subset of cells where we are performing the analysis, i.e. E1, E1col1, or full')
parser.add_argument('tkey', type=int, default=0, help='whether to use days as `tkey` when computing velocity graph')
args = parser.parse_args()

if not args.tkey:
    data_path = 'data/{}_velocity.h5ad'.format(args.subset)
else:
    data_path = 'data/{}_velocity_tkey.h5ad'.format(args.subset)

adata = sc.read(data_path)

scv.pl.velocity_embedding_stream(adata, basis='umap', legend_fontsize=12, dpi=200, title='{} by type'.format(args.subset), color='type', 
    save='{}_stream_type_tkey{}.png'.format(args.subset, args.tkey))

scv.pl.velocity_embedding_stream(adata, basis='umap', legend_fontsize=12, dpi=200, title='{} by day'.format(args.subset), color='day', 
    save='{}_stream_day_tkey{}.png'.format(args.subset, args.tkey))

scv.pl.velocity_embedding(adata, basis='umap', dpi=200, density=0.75, color='type', 
    save='{}_arrow_type_tkey{}.png'.format(args.subset, args.tkey))

scv.pl.velocity_embedding(adata, basis='umap', dpi=200, density=0.75, color='day', 
    save='{}_arrow_day_tkey{}.png'.format(args.subset, args.tkey))