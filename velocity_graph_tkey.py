import scanpy as sc
import scvelo as scv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('subset', help='the subset of cells where we are performing the analysis, i.e. E1, E1col1, or full')
args = parser.parse_args()

data_path = 'data/{}_velocity.h5ad'.format(args.subset)

adata = sc.read(data_path)
scv.tl.velocity_graph(adata, tkey='day')
adata.write('data/{}_velocity_tkey.h5ad'.format(args.subset))