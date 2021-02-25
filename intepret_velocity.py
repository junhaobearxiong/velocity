import scanpy as sc
import scvelo as scv
import argparse
import numpy as np

np.random.seed(1)

parser = argparse.ArgumentParser()
parser.add_argument('subset', help='the subset of cells where we are performing the analysis, i.e. E1, E1col1, or full')
parser.add_argument('plot', help='which type of plot to produce')
parser.add_argument('given_root', type=int, default=0, help='whether root cells are given (relevant for latent time or velocity graph')
parser.add_argument('tkey', type=int, default=0, help='whether to use days as `tkey` when computing velocity graph')
args = parser.parse_args()

if not args.tkey:
    data_path = 'data/{}_velocity.h5ad'.format(args.subset)
else:
    data_path = 'data/{}_velocity_tkey.h5ad'.format(args.subset)

adata = sc.read(data_path)

if args.plot == 'latent_time':
    top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
    latent_time_save_path = '{}_latent_time_root{}_tkey{}.png'.format(args.subset, args.given_root, args.tkey)
    top_genes_save_path = '{}_top_genes_by_latent_time_root{}_tkey{}.png'.format(args.subset, args.given_root, args.tkey)
    if not args.given_root:
        # estimate gene-shared latent time
        scv.tl.latent_time(adata)
    else:
        # root cells given by day 0 cells
        adata.obs['is_day0'] = (adata.obs['day'] == 0).astype(int)
        scv.tl.latent_time(adata, root_key='is_day0')

    # plot embedding color by latent time
    scv.pl.scatter(adata, color='latent_time', cmap='gnuplot', dpi=200,
        save=latent_time_save_path)
    # plot top genes expression ordered by pseudotime
    scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='type', n_convolve=100, 
        save=top_genes_save_path)


elif args.plot == 'velocity_graph':
    scv.pl.velocity_graph(adata, color='type', dpi=200, n_neighbors=5,
        save='{}_{}_type_tkey{}.png'.format(args.subset, args.plot, args.tkey))
    scv.pl.velocity_graph(adata, color='day', dpi=200, n_neighbors=5,
        save='{}_{}_day_tkey{}.png'.format(args.subset, args.plot, args.tkey))


elif args.plot == 'velocity_pseudotime':
    velocity_pseudotime_save_path = '{}_velocity_pseudotime_root{}_tkey{}.png'.format(args.subset, args.given_root, args.tkey)
    paga_save_path = '{}_paga_root{}_tkey{}.png'.format(args.subset, args.given_root, args.tkey)
    if not args.given_root:
        scv.tl.velocity_pseudotime(adata)
    else:
        # randomly select 1 day 0 ipsc as root cell
        root_indices = np.arange(adata.shape[0])[(adata.obs['day'] == 0) & (adata.obs['type'] == 'IPSC')]
        root_idx = np.random.choice(root_indices)
        scv.tl.velocity_pseudotime(adata, root_key=root_idx)

    scv.pl.scatter(adata, color='velocity_pseudotime', color_map='gnuplot', dpi=200, 
        save=velocity_pseudotime_save_path)

    # this is needed due to a current bug - bugfix is coming soon.
    adata.uns['neighbors']['distances'] = adata.obsp['distances']
    adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

    # plot paga graph extended by velocity-inferred directionality
    # using `velocity_pseudotime` as prior
    scv.tl.paga(adata, groups='type')
    df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
    df.style.background_gradient(cmap='Blues').format('{:.2g}')

    scv.pl.paga(adata, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5, dpi=150,
        save=paga_save_path)


elif args.plot == 'top_genes':
    # plot top genes phase portrait
    top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index
    scv.pl.scatter(adata, basis=top_genes[:15], ncols=5, frameon=False, dpi=150, color='type',
        save='{}_top_genes_phase_portrait.png'.format(args.subset))


elif args.plot == 'speed_coherence':
    scv.tl.velocity_confidence(adata)
    keys = 'velocity_length', 'velocity_confidence'
    scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], dpi=200,
        save='{}_{}.png'.format(args.subset, args.plot))
    df = adata.obs.groupby('type')[keys].mean().T
    df.style.background_gradient(cmap='coolwarm', axis=1)


else:
    raise ValueError('plot {} is not implemented'.format(args.plot))
