import os
import re
import os.path
import pandas as pd
import numpy as np
import scanpy as sc
import scanpy.external as sce
from scipy.sparse import csr_matrix, isspmatrix, isspmatrix_csr, issparse


def del_raw(adata):
    if adata.raw:
        del(adata.raw)

def enforce_sparsity(adata):
    if isspmatrix(adata.X):
        if not isspmatrix_csr(adata.X):
            adata.X = adata.X.tocsr()
    else:
        adata.X = csr_matrix(adata.X)

def fix_nan_x(adata):
    if np.isnan(adata.X.data).any():
        adata.X.data = np.select(
                [np.isnan(adata.X.data), np.full_like(adata.X.data, True, dtype=bool)],
                [np.full_like(adata.X.data, 0.0), adata.X.data])
        
def filter_nan_var_names(adata):
    return adata[:, adata.var_names.notnull()].copy()

def ensure_unique_idx(adata):
    for i in ["obs_names", "var_names"]:
        if not getattr(adata, i).is_unique:
            getattr(adata, f"{i}_make_unique")()

def calc_n_genes(adata):
    if "n_genes" not in adata.obs:
        if issparse(adata.X):
            adata.obs["n_genes"] = np.sum(adata.X > 0, axis=1).A1
        else:
            adata.obs["n_genes"] = np.sum(adata.X > 0, axis=1)

def find_doublets(adata, **kwargs):
    sc.pp.scrublet(adata, **kwargs)
    doublet_counts = adata.obs["predicted_doublet"].value_counts()
    
    return adata, doublet_counts

def filter_low_quality_cells(adata, min_genes=None, min_cells=None, min_counts=None, max_counts=None, 
                             min_mt_fraction=None, max_mt_fraction=None, use_qc_metrics=True):
    if use_qc_metrics:
        adata.var['mt'] = adata.var_names.str.contains(r'^(?i)mt-', regex=True)
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

    if min_genes:
        sc.pp.filter_cells(adata, min_genes=min_genes)
    
    if min_cells:
        sc.pp.filter_genes(adata, min_cells=min_cells)

    if min_counts or max_counts:
        if min_counts:
            sc.pp.filter_cells(adata, min_counts=min_counts)
        if max_counts:
            adata = adata[adata.obs['total_counts'] <= max_counts]

    if use_qc_metrics and min_mt_fraction is not None:
        adata = adata[adata.obs['pct_counts_mt'] >= min_mt_fraction]

    if use_qc_metrics and max_mt_fraction is not None:
        adata = adata[adata.obs['pct_counts_mt'] <= max_mt_fraction]

def normalize_total(adata, target_sum=1e4):
    if "n_counts" not in adata.obs:
        adata.obs["n_counts"] = np.ravel(np.sum(adata.X, axis=1))
    sc.pp.normalize_total(adata, target_sum=target_sum)

def standardize_gene_symbols(adata, species):
    species = species.capitalize()
    
    transformations = {
        "Human": lambda x: x.upper(),
        "Mouse": lambda x: x.capitalize()
    }

    if species not in transformations:
        raise ValueError(f"Species '{species}' is not supported. Available options: {list(transformations.keys())}")

    adata.var_names = [transformations[species](x) for x in adata.var_names]
    