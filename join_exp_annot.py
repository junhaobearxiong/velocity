'''
Join the velocyto loom file with the expression annotation
and find the set of cells with both spliced/unspliced read counts
and cell type/day annotation
'''

import scvelo as scv
import scanpy as sc
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('subset', help='the subset of cells where we are performing the analysis, i.e. E1, E1col1, or full')
args = parser.parse_args()

input_exp_path = 'data/seurat.annotated.sct.h5ad'
input_velo_path = 'data/{}.loom'.format(args.subset)
output_path = 'data/{}_joint.h5ad'.format(args.subset)

adata = scv.read(input_velo_path)
scv.utils.show_proportions(adata)
print('velocyto data shape: {}'.format(adata.shape))
adata_exp = sc.read(input_exp_path)
print('expression data shape: {}'.format(adata_exp.shape))

def process_velocyto_cellname(x):
    strings = x.split('_')
    return strings[1] + '_' + strings[3].split(':')[1]

def process_expression_cellname(x):
    strings = x.split('_')
    return strings[1] + '_' + strings[2]

new_velo_index = pd.Series(adata.obs.index).apply(lambda x: process_velocyto_cellname(x))
new_exp_index = pd.Series(adata_exp.obs.index).apply(lambda x: process_expression_cellname(x))
adata.obs.index = new_velo_index
adata_exp.obs.index = new_exp_index

# find the cells present in both velocity file and expression file
joint_index = adata_exp.obs.join(adata.obs, how='inner').index
adata = adata[joint_index]
adata_exp = adata_exp[joint_index]

# combine the annotation of expression data with the velocity data
adata.obs = adata.obs.join(adata_exp.obs)
adata.obs['day'] = adata.obs['diffday'].apply(lambda x: int(x[3:]))

adata.write(output_path)