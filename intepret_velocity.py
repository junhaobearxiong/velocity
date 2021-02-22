import scanpy as sc
import scvelo as scv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('subset', help='the subset of cells where we are performing the analysis, i.e. E1, E1col1, or full')
parser.add_argument('plot', help='which type of plot to produce')
args = parser.parse_args()

data_path = 'data/{}_velocity.h5ad'.format(args.subset)

adata = sc.read(data_path)

if args.plot == 'latent_time':
    # estimate gene-shared latent time
    scv.tl.latent_time(adata)
    # plot embedding color by latent time
    scv.pl.scatter(adata, color='latent_time', cmap='gnuplot', dpi=200,
        save='{}_latent_time.png'.format(args.subset))
    # plot top genes expression ordered by pseudotime
    top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
    scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='type', n_convolve=100, 
        save='{}_top_genes_by_latent_time.png'.format(args.subset))

elif args.plot == 'latent_time_with_root':
    adata.obs['is_day0'] = (adata.obs['day'] == 0).astype(int)
    # estimate gene-shared latent time with root cells given by day 0 cells
    scv.tl.latent_time(adata, root_key='is_day0')
    # plot embedding color by latent time
    scv.pl.scatter(adata, color='latent_time', cmap='gnuplot', dpi=200,
        save='{}_latent_time_with root.png'.format(args.subset))
    # plot top genes expression ordered by pseudotime
    top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
    scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='type', n_convolve=100, 
        save='{}_top_genes_by_latent_time_with_root.png'.format(args.subset))

elif args.plot == 'top_genes':
    # plot top genes phase portrait
    top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index
    scv.pl.scatter(adata, basis=top_genes[:15], ncols=5, frameon=False, dpi=150, color='type',
        save='{}_top_genes_phase_portrait.png'.format(args.subset))

elif args.plot == 'velocity_graph':
    scv.tl.velocity_graph(adata, tkey='day')
    scv.pl.velocity_graph(adata, color='type', threshold=.1, dpi=200, n_neighbors=5,
        save='{}_{}_type.png'.format(args.subset, args.plot))
    scv.pl.velocity_graph(adata, color='day', threshold=.1, dpi=200, n_neighbors=5,
        save='{}_{}_day.png'.format(args.subset, args.plot))

    # plot embedding colored by velocity pseudotime, which is based on velocity graph
    scv.tl.velocity_pseudotime(adata)
    scv.pl.scatter(adata, color='velocity_pseudotime', color_map='gnuplot', dpi=200, 
        save='{}_velocity_pseudotime.png'.format(args.subset))

    # this is needed due to a current bug - bugfix is coming soon.
    adata.uns['neighbors']['distances'] = adata.obsp['distances']
    adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

    # plot paga graph extended by velocity-inferred directionality
    scv.tl.paga(adata, groups='type')
    df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
    df.style.background_gradient(cmap='Blues').format('{:.2g}')

    scv.pl.paga(adata, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5, dpi=150,
        save='{}_paga.png'.format(args.subset))

elif args.plot == 'speed_coherence':
    scv.tl.velocity_confidence(adata)
    keys = 'velocity_length', 'velocity_confidence'
    scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], dpi=200,
        save='{}_{}.png'.format(args.subset, args.plot))
    df = adata.obs.groupby('type')[keys].mean().T
    df.style.background_gradient(cmap='coolwarm', axis=1)

else:
    raise ValueError('plot {} is not implemented'.format(args.plot))
