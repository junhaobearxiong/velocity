import scanpy as sc
import scvelo as scv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('subset', help='the subset of cells where we are performing the analysis, i.e. E1, E1col1, or full')
args = parser.parse_args()

data_path = 'data/{}_velocity.h5ad'.format(args.subset)

adata = sc.read(data_path)
scv.pl.velocity_embedding_stream(adata, basis='umap', legend_fontsize=12, title='', smooth=.8, min_mass=4, color='type', save='stream_type.png')

scv.pl.velocity_embedding_stream(adata, basis='umap', legend_fontsize=12, title='', smooth=.8, min_mass=4, color='day', save='stream_day.png')

scv.pl.velocity_embedding(adata, basis='umap', arrow_length=2, arrow_size=2, dpi=200, color='type', save='arrow_type.png')

scv.pl.velocity_embedding(adata, basis='umap', arrow_length=2, arrow_size=2, dpi=200, color='day', save='arrow_day.png')