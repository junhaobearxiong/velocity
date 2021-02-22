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
    # estimate pseudotime
    scv.tl.latent_time(adata)
    # plot embedding color by pseudotime
    scv.pl.scatter(adata, color='latent_time', cmap='gnuplot', dpi=200,
        save='{}_latent_time.png'.format(args.subset))
    # plot top genes expression ordered by pseudotime
    top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
    scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='type', n_convolve=100, 
        save='{}_top_genes_by_time.png'.format(args.subset))

elif args.plot == 'top_genes':
    # plot top genes phase portrait
    top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index
    scv.pl.scatter(adata, basis=top_genes[:15], ncols=5, frameon=False, dpi=150, color='type',
        save='{}_top_genes_phase_portrait.png'.format(args.subset))

elif args.plot == 'velocity_graph':
    scv.pl.velocity_graph(adata, color='type', dpi=200, threshold=.1, save='{}_{}.png'.format(args.subset,))

elif args.plot == 'speed_coherence':
    scv.tl.velocity_confidence(adata)
    keys = 'velocity_length', 'velocity_confidence'
    scv.pl.scatter(adata, c=keys, cmap='coolwarm', dpi=200, perc=[5, 95], save='{}_{}.png'.format(args.subset, args.plot))
    df = adata.obs.groupby('type')[keys].mean().T
    df.style.background_gradient(cmap='coolwarm', axis=1)

else:
    raise ValueError('plot {} is not implemented'.format(args.plot))
