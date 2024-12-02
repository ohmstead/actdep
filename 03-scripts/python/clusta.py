"""
This script has some tools that are useful for pre-processing and clustering my
snRNA-seq data.

In general, the workflow relies heavily on steps taken from Parse's scanpy
tutorial, which can be viewed here (requires Parse account login to their support suite):
https://support.parsebiosciences.com/hc/en-us/articles/360052794312-Scanpy-Tutorial-65k-PBMCs.
"""

import os
import warnings
import numpy as np
import scipy as sp
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.io as sio
import matplotlib.pyplot as plt
from scranPY import compute_sum_factors
from datetime import datetime

warnings.filterwarnings("ignore")

plt.rcParams['font.family'] = 'avenir'

sc.settings.set_figure_params(dpi=150, fontsize=10, dpi_save=300, figsize=(5,5), format='png')


def get_data_path(sublibrary, dataset):
    """
    Returns path to folder containing sublibrary data. Used for reading/writing AnnData objects.

    Parameters
    ----------
    sublibrary : string
        MUST be "sublibrary1", "sublibrary2", or "combined"

    Returns
    -------
    data_path : string
        string for path leading to folder containing sublibrary data
    """

    # validate input
    if sublibrary not in ['sublibrary1', 'sublibrary2', 'combined']:
        raise ValueError('sublibrary must be "sublibrary1", "sublibrary2", or "combined"')
    
    # return path based on location of seq directory
    if dataset == 'Dec2023':
        parent_location = '/Volumes/jack/seq/analysis_Dec2023/'
    elif dataset == 'Feb2024':
        parent_location = '/Volumes/jack/seq/analysis_Feb2024/'
    elif dataset == 'May2024':
        parent_location = '/Volumes/jack/seq/analysis_May2024/'
    elif dataset == 'merge':
        parent_location = '/Volumes/jack/seq/analysis_merge/'
    else:
        raise ValueError('dataset in clusta.get_data_path() must be "Dec2023", "Feb2024", or "merge"')
    
    data_path = os.path.join(parent_location, sublibrary, '0_all-sample/DGE_filtered/')

    return data_path


def read_raw_count_matrix(sublibrary, dataset):
    """
    Generates an AnnData object from raw .mtx file output by Parse pipeline.
    Populates adata.obs with cell metadata and adata.var with gene metadata.

    Relies on outputs being formatted according to outputs from the Parse
    pipeline.

    Parameters
    ----------
    sublibrary : string
        MUST be "sublibrary1", "sublibrary2", or "combined"

    Returns
    -------
    adata : AnnData object
        AnnData object containing raw count matrix and cell/gene metadata
    """

    data_path = get_data_path(sublibrary, dataset)

    # validate sublibrary input
    if sublibrary not in ['sublibrary1', 'sublibrary2', 'combined']:
        raise ValueError('sublibrary must be "sublibrary1", "sublibrary2", or "combined"')

    # import data
    adata = sc.read_mtx(data_path + 'count_matrix.mtx')

    # read in gene and cell data
    gene_data = pd.read_csv(data_path + 'all_genes.csv')
    cell_meta = pd.read_csv(data_path + 'cell_metadata.csv')

    # find genes with nan values and filter
    gene_data = gene_data[gene_data.gene_name.notnull()]
    notNa = gene_data.index
    notNa = notNa.to_list()

    # remove genes with nan values and assign gene names
    adata = adata[:,notNa]
    adata.var = gene_data
    adata.var.set_index('gene_name', inplace=True)
    adata.var.index.name = None
    adata.var_names_make_unique()

    # add cell meta data to anndata object
    adata.obs = cell_meta
    adata.obs.set_index('bc_wells', inplace=True)
    adata.obs.index.name = None
    adata.obs_names_make_unique()

    # add a sublibrary column
    if sublibrary == 'combined':
        adata.obs['sublibrary'] = np.where(adata.obs.index.str.endswith('s1'), 'sublibrary1', 'sublibrary2')
        adata.uns['sublibrary'] = 'combined'
    else:
        adata.obs['sublibrary'] = sublibrary
        adata.uns['sublibrary'] = sublibrary
    
    adata.uns['dataset'] = dataset

    return adata


def qc_filter_cells(adata, plots_on=True, 
                    min_counts=1500,
                    max_counts=80000, 
                    min_genes=1500, 
                    max_genes=8000, 
                    max_mito=1
                    ):
    """
    Filters cells based on transcript count, gene count, and mitochondrial
    reads. Plots QC metrics if desired.

    Parameters
    ----------
    adata : AnnData object with raw count matrix
        has raw counts
    plots_on : bool, optional
        show plots, by default True
    min_counts : int, optional
        minimum number of unique transcripts for a cell to be included, by
        default 1500
    max_counts : int, optional
        maximum number of unique transcripts for a cell to be included, by
        default 80000
    min_genes : int, optional
        mininum number of unique genes for a cell to be included, by default
        1500
    max_genes : int, optional
        maximum number of unique genes for a cell to be included, by default
        8000
    max_mito : int, optional
        maximum percent allowed of mitochrondrial reads in one cell, by default
        1

    Returns
    -------
    adata : AnnData object
        contains only filtered cells
    """

    # calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    if plots_on:
        sc.pl.violin(adata, ['n_genes_by_counts'], save='_n_genes', jitter=0.4)
        sc.pl.violin(adata, ['total_counts'], save='_total_counts', jitter=0.4, log=True)
        sc.pl.violin(adata, ['pct_counts_mt'], save='_mito_pct', jitter=0.4)

    # exclude cells based on scalar thresholds
    sc.pp.filter_genes(adata, min_cells=5)  # remove genes expressed in fewer than 5 cells
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_cells(adata, max_counts=max_counts)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, max_genes=max_genes)
    adata = adata[adata.obs.pct_counts_mt < max_mito,:]

    if plots_on:
        sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save='_gene_vs_transcript_counts')

    print('median transcript count per cell: ' + str(adata.obs['tscp_count'].median(0)))
    print('median gene count per cell: ' + str(adata.obs['gene_count'].median(0)))
    print(adata.shape) # print number of cells remaining

    return adata


def normalize_and_log1p_transform(adata):
    """
    Normalizes data using scran's computeSumFactors function. Afterward, counts
    are log1p-transformed.
    
    Raw counts are stored in adata.layers['rawcounts'].

    Parameters
    ----------
    adata : AnnData object
        contains raw counts

    Returns
    -------
    adata : AnnData object
        contains normalized, log1p-transformed counts
    """

   # save raw counts before normalizing
    adata.layers['rawcounts'] = adata.X

    # convert sparse to dense matrix; scran only works on dense
    adata.layers['denoised'] = adata.X.A
    adata.X = adata.layers['denoised']
    sp.sparse.issparse(adata.X)

    # scran normalization 
    compute_sum_factors(adata=adata, normalize_counts=True, log1p=True)

    adata.X = sp.sparse.csr_matrix(adata.layers['scranPY']).copy()
    adata.X.shape

    return adata


def select_hvgs(adata,
                save_new_adata=False,
                dispersion_min_mean=0.0125,
                dispersion_max_mean=3,
                dispersion_min=0.25,):
    """
    Selects highly variable genes based on dispersion. All other genes are
    removed from AnnData object.

    scanpy documentation:
    https://scanpy.readthedocs.io/en/stable/api/scanpy.pp.highly_variable_genes.html

    Parameters
    ----------
    adata : AnnData object
        AnnData with normalized, log1p-transformed counts in adata.X
    dispersion_min_mean : float, optional
        see scanpy documentation, by default 0.0125
    dispersion_max_mean : int, optional
        see scanpy documentation, by default 3
    dispersion_min : float, optional
        see scanpy documentation, by default 0.25

    Returns
    -------
    adata
        AnnData object with only highly variable genes
    """

    # get data path
    data_path = get_data_path(adata.uns['sublibrary'] , adata.uns['dataset'])

    # find highly variable genes based on dispersion
    if adata.uns['sublibrary'] == 'combined':
        batch = 'sublibrary'
    else:
        batch = None
    sc.pp.highly_variable_genes(adata,
                                min_mean=dispersion_min_mean, 
                                max_mean=dispersion_max_mean, 
                                min_disp=dispersion_min, 
                                batch_key=batch)
    sc.pl.highly_variable_genes(adata, save='')

    # save raw expression values before variable gene subset
    if save_new_adata:
        adata.raw = adata
        adata.write(os.path.join(data_path, 'adata_preFeatSelection.h5ad'))

    # remove low-variance genes
    adata = adata[:, adata.var.highly_variable]

    return adata


def regress_out_covariates(adata, covariates=['total_counts']):
    """
    Regresses out covariates from data. Covariates are specified in the
    function. 

    Parameters
    ----------
    adata : AnnData object
        AnnData object with normalized, log1p-transformed counts in adata.X

    Returns
    -------
    adata : AnnData object
        AnnData object with covariates regressed out
    """

    # validate that covariates are a list of strings
    if not isinstance(covariates, list):
        raise TypeError('covariates must be a list of strings')

    # regress out covariates
    sc.pp.regress_out(adata, keys=covariates)

    return adata


def scale_data(adata):
    """
    Scales data to unit variance and zero mean.

    Parameters
    ----------
    adata : AnnData object
        AnnData object with normalized, log1p-transformed counts in adata.X

    Returns
    -------
    adata : AnnData object
        AnnData object with scaled data in adata.X
    """

    sc.pp.scale(adata)

    return adata


def run_pca(adata, n_pcs=50):
    """
    Runs PCA on data. Saves PCA coordinates in adata.obsm['X_pca'].

    Parameters
    ----------
    adata : AnnData object
        AnnData object with scaled data in adata.X
    n_pcs : int, optional
        number of principal components to compute, by default 50

    Returns
    -------
    adata : AnnData object
        AnnData object with PCA coordinates in adata.obsm['X_pca']
    """

    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_pcs)
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=n_pcs, save='')

    return adata


def make_neighborhood_graph(adata, n_neighbors=15, num_PCs_to_use=50):
    """
    Makes a neighborhood graph for clustering.

    If both sublibraries are included, it uses the BBKNN method to perform batch
    correction (see below for documentation):
    https://scanpy.readthedocs.io/en/stable/generated/scanpy.external.pp.bbknn.html

    Parameters
    ----------
    adata : AnnData object
        AnnData object with PCA coordinates in adata.obsm['X_pca']
    n_neighbors : int, optional
        number of neighbors to use in graph, by default 15
    num_PCs_to_use : int, optional
        number of principal components to use in graph, by default 50

    Returns
    -------
    adata : AnnData object
        AnnData object with neighborhood graph in adata.uns['neighbors']
    """

    if adata.uns['sublibrary'] == 'combined':
        sc.external.pp.bbknn(adata, batch_key='sublibrary')
    else:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=num_PCs_to_use)

    return adata


def louvain_cluster(adata, resolution=0.5):
    """
    Performs Louvain clustering on neighborhood graph.

    Parameters
    ----------
    adata : AnnData object
        AnnData object with neighborhood graph in adata.uns['neighbors']
    resolution : float, optional
        resolution parameter for Louvain clustering, by default 1.0

    Returns
    -------
    adata : AnnData object
        AnnData object with Louvain clusters in adata.obs['louvain']
    """

    sc.tl.umap(adata)
    sc.tl.louvain(adata, resolution=resolution)

    sc.pl.umap(adata, color='louvain')
    sc.pl.umap(adata, color='sublibrary')

    return adata


def make_adata_lean(adata, adata_filename):
    """
    Removes unnecessary data from AnnData object to make it leaner for use with
    Allen Brain Institute's MapMyCells tool. Writes leaner AnnData object to file.

    MapMyCells requires raw counts, so adata.X is replaced with adata.layers['rawcounts'].

    Parameters
    ----------
    adata : AnnData object
        AnnData object that has been filtered and clustered
    adata_filename : string
        the name of the file to save the leaner AnnData object as
    """

    # check adata_filename ends with .h5ad
    if not adata_filename.endswith('.h5ad'):
        raise ValueError('adata_filename must end with ".h5ad"')
    
    data_path = get_data_path(adata.uns['sublibrary'], adata.uns['dataset'])
    adata_lean_path = os.path.join(data_path, adata_filename)

    # check if adata_filename already exists and ask user to confirm overwrite
    if os.path.exists(adata_lean_path):
        overwrite = input(f"The file {adata_filename} already exists. Do you want to overwrite it? (y/n): ")
        if overwrite.lower() != 'y':
            print("File not overwritten. Exiting...")
            return

    # trim to bare-bones metadata and save as 'adata_lean' for use w/ MapMyCells
    adata_lean = ad.AnnData(adata.layers['rawcounts'])

    adata_lean.var_names = adata.var_names
    adata_lean.obs_names = adata.obs_names

    # lastly, store sublib info
    adata_lean.uns['sublibrary'] = adata.uns['sublibrary']

    # write lean data to file
    adata_lean.write(adata_lean_path, compression='gzip')


def merge_allen_celltypes(adata, allen_csv_path):
    """
    Merges Allen Brain Institute cell types with adata.obs.

    Parameters
    ----------
    adata : AnnData object
        AnnData object with Louvain clusters in adata.obs['louvain']
    allen_csv_path : string
        path to csv file containing Allen Brain Institute cell types

    Returns
    -------
    adata : AnnData object
        AnnData object with Allen Brain Institute cell types in adata.obs['allen_celltypes']
    """

    # confirm that allen_csv_path exists in the correct sublibrary folder
    data_path = get_data_path(adata.uns['sublibrary'], adata.uns['dataset'])
    if not os.path.exists(os.path.join(data_path, allen_csv_path)):
        raise ValueError(f"allen_csv_path does not exist in {data_path}")

    # read in metadata from two sources: (1) MapMyCells (2) from adata object produced by preprocessing
    meta = pd.read_csv(os.path.join(data_path, allen_csv_path), header=4)
    meta.set_index(keys='cell_id', drop=True, inplace=True)
    meta.index.name = None

    # merge MapMyCells metadata into existing adata object
    adata.obs = adata.obs.join(meta, how='inner', validate='1:1')

    return adata


def get_taxonomy_cmap(adata):
    """
    Loads taxonomy information from the Allen Brain Atlas and creates a
    consistent color map for each taxon level. This means that class, subclass,
    supertype, and clusters will all have consistent coloring across AnnData
    objects and clustering parameters.

    Parameters
    ----------
    adata : AnnData object
        should have Allen cell type assignment in adata.obs

    Returns
    -------
    adata: AnnData object
        adata now has color maps for each taxon level in
        adata.uns['class_name_colors'], etc.
    """

    from collections import OrderedDict

    # check that Allen cell types have been mapped to adata
    if 'class_name' not in adata.obs.columns:
        raise ValueError('adata.obs must have Allen cell types in adata.obs["class_name"]')
    
    # load in information about taxonomy from Allen Brain Atlas
    allen_taxonomy = pd.read_excel('/Users/jack/R/seq/cl.df_CCN202307220.xlsx')

    # for each taxon level, make a dictionary mapping cell types to colors
    taxa = ['class', 'subclass', 'supertype', 'cluster']
    
    for taxon in taxa:
        adata_key = taxon + '_name'
        
        # create dict for all possible types in the Allen nomenclature
        taxon_names = allen_taxonomy[taxon + '_id_label'].unique()
        taxon_colors = dict(zip(taxon_names, gb.create_palette(len(taxon_names))))

        # throw out cell types that don't occur actually appear in adata
        taxon_colors = {k: v for k, v in taxon_colors.items() if k in adata.obs[adata_key].astype('category').cat.categories}

        # order the dict so that the colors are in the same order as a sorted list of cell types
        taxon_colors_ordered = dict((k, taxon_colors[k]) for k in adata.obs[adata_key].astype('category').cat.categories)
        adata.uns[taxon + '_colorDict'] = taxon_colors_ordered
        
        taxon_colors_ordered = list(taxon_colors_ordered.values())
        adata.uns[adata_key + '_colors'] = taxon_colors_ordered

    return adata


def merge_small_clusters(adata, min_cluster_size=2, taxa=['class', 'subclass', 'supertype', 'cluster']):
    """
    Merge small clusters in the given AnnData object.

    Parameters
    ----------
    adata : AnnData object
        The AnnData object containing the taxa to be merged.
    min_cluster_size : int, optional
        The minimum number of cells required for a cluster to be considered as a
        small cluster. Default is 2.
    taxa : list, optional
        The taxonomic levels to consider for merging clusters. Default is for
        all Allen taxonomic levels: ['class', 'subclass', 'supertype', 'cluster'].

    Returns
    -------
    adata : AnnData object
        The AnnData object with merged small clusters.

    Notes
    -----
    This function re-assigns all clusters with a number of cells less than
    `min_cluster_size` to a cluster named 'SMALL-x', where 'x' is the taxonomic
    level. Importantly, this retains its bootstrapping probability to the
    original taxon.
    """
    
    # re-assign all classes with 1 cell to class 'small'
    for taxon in taxa:
        clusts = adata.obs[taxon + '_name'].value_counts()

        # replace taxon name/label to 'SMALL-x', where x is the taxon level
        for clust in clusts.index:
            if clusts[clust] < min_cluster_size:
                adata.obs[taxon + '_name'].replace(clust, 'SMALL-' + taxon, inplace=True)
                adata.obs[taxon + '_label'].replace(clust, 'SMALL-' + taxon, inplace=True)
    
    return adata