"""
The whole purpose of this script is to make the DGE analysis I did using MAST reproducible in R.

To do this, I want to make a script that has the python code I used to set up the adata object that eventually was input into MAST.

This is to ensure that I'm running MAST on identical cells and genes as I did in python. Just in R.

The code I'm drawing from is from dge.ipynb. I'm going to copy the code from that notebook into this script.
"""

import warnings

warnings.filterwarnings("ignore")

import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import clusta
import utils
import os

sc.settings.verbosity = 0

def get_adata(DGE_folder):
    # import adata
    adata = ad.read_h5ad(os.path.join(DGE_folder, 'adata.h5ad'))

    # cut down to just subclass 016 CA1-ProS Glut
    adata = adata[adata.obs['subclass_name'] == '037 DG Glut']

    # get cell IDs from obs
    ca1_ids = adata.obs.index.tolist()

    # detect 'combined', 'sublibrary1', 'sublibrary2' in DGE_folder string
    if 'combined' in DGE_folder:
        dset = 'combined'
    elif 'sublibrary1' in DGE_folder:
        dset = 'sublibrary1'
    elif 'sublibrary2' in DGE_folder:
        dset = 'sublibrary2'


    # load in raw data from cell-count mat
    adata_ca1 = clusta.read_raw_count_matrix(dset)
    adata_ca1

    # only keep cells in cell_ids
    adata_ca1 = adata_ca1[adata_ca1.obs.index.isin(ca1_ids)]
    adata_ca1.obs['condition'] = adata_ca1.obs['sample'].str[:2]

    adata_ca1

    # sc.pp.normalize_total(adata_ca1, target_sum=1e4)
    # sc.pp.log1p(adata_ca1)
    # sc.pp.filter_genes(adata_ca1, min_cells=30)

    adata_ca1 = utils.prep_anndata(adata_ca1)

    # make adata_ca1.obs['condition'] the first 2 letters of the 'sample' column
    adata_ca1.obs = pd.concat([adata_ca1.obs, adata.obs], axis=1)
    # adata_ca1.obs['subclass_name'] = 'CA1'

    return adata_ca1