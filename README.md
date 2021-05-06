# RNA velocity on (time-course) single-cell expression data

## Intro
See the [slides](https://docs.google.com/presentation/d/1_Rd5fgE_AYLVqPdnvD-ux-t5Ey3bKZNE5_dChNtNkTk/edit?usp=sharing) here for the workflow and results of an example analysis.

## Preparation 
First, follow documentation on `scanpy` ([link](https://scanpy.readthedocs.io/en/stable/installation.html)) and `scvelo` ([link](https://scvelo.readthedocs.io/installation.html)) to install the packages and their dependencites. You can also install `cellrank` ([link](https://cellrank.readthedocs.io/en/stable/installation.html)), but it's not necessary for running most of the scripts. 

This repository was built with the following version of the main packages:
```
cellrank==1.2.0
scanpy==1.7.0
scvelo==0.2.3
```

Then, create the following directories under the main directory: `data/` (where the main input data should be saved, and the preprocessed data will also be saved here) and `figures/` (where the figures will be saved).

## Workflow
Below is a standard workflow for preprocessing, estimating and analyzing RNA velocity on (time-course) single-cell expression data. It is largely based on the [tutorial](https://scvelo.readthedocs.io/getting_started.html) on the `scvelo` website.

#### Run `velocyto`
Starting with the `bam` files, we first need to run the command line interface (CLI) of `velocyto` to generate unspliced/spliced count matrices (see [tutorial](http://velocyto.org/velocyto.py/tutorial/cli.html) for many more details). For the example analysis, I first used `samtools` to sort each bam files by cell barcodes (as required by `velocyto`), then use the `run` command of the CLI. Users can modify `sort_bam.sh` and `run_velocyto.sh` for this step. Note that you may need to create additional directories, change file names, download additional files, etc. 

#### Preprocess data before running `scvelo`
For the example analysis, since I have a `loom` file for each round, collection, collection day, I need to combine them first. This is easily done with the `loompy` package. Users can modify `combine_loom.py` to combine all the relevant `loom` files into a single `loom` file.

Then, I filtered out certain cells by taking the intersection between the `loom` file and an already-filtered, annotated expression file. This also adds annotations of cell types and time points to the `loom` files. I also filtered out genes, normalized the data, computed neighbours and moments, estimated PCA and UMAP, etc. Users can modify `filter_cells_genes.py` for this step.

#### Estimate velocity with `scvelo`
Taking the output from `filter_cells_genes.py`, users can run `estimate_velocity.py` to estimate RNA velocity.

#### Analyze velocity 
It's often useful to project and visualize velocity on a lower-dimensional embedding. `project_velocity.py` generates such figures.

There are many other analyses one can perform with the velocity estimates using `scvelo`. Some of these visualizations can be done with `intepret_velocity.py` by feeding in the appropriate command line arguments. 

#### Additional analyses 
When estimating the velocity graph, `scvelo` allows an additional argument `tkey` for time-course data. It's unclear if it does the right thing (see github [issue](https://github.com/theislab/scvelo/issues/367)), but users can get velocity estimates under the `tkey` prior using `velocity_graph_tkey.py`.

In the example analysis we have multiple samples / cell lines. After running the analysis on all the data, I split the `loom` files by cell lines and generate plots for each cell line using `create_file_by_id.py` and `project_velocity_by_line.py`.

#### Run `cellrank`
I did not have success getting intepretable results using `cellrank`, but users can use `run_cellrank.py` to start the analysis with `cellrank`. 
