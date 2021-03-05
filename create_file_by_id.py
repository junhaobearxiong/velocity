'''
Create a separate h5ad file for each cell line
'''

import scvelo as scv
import argparse

input_path = 'data/full_joint.h5ad'
adata = scv.read(input_path, cache=True)

for ind in adata.obs['individual'].unique():
    output_adata = adata[adata.obs['individual'] == ind]
    output_path = 'data/{}_joint.h5ad'.format(ind)
    output_adata.write(output_path)