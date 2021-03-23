import scvelo as scv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('data1')
parser.add_argument('data2')
args = parser.parse_args()

adata1 = scv.read(args.data1, cache=True)
adata2 = scv.read(args.data2, cache=True)

genes1 = adata1.var.index
genes2 = adata2.var.index
print('# genes in data 1: {}'.format(genes1.size))
print('# genes in data 2: {}'.format(genes2.size))
print('# genes in intersection: {}'.format(genes1.intersection(genes2).size))