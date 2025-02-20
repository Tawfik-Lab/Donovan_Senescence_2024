import scanpy as sc
import numpy as np
import pandas as pd
import scipy.stats as sps
from scipy.sparse import issparse
from util import adata_filter, adata_filter_mask, get_geneset
from statsmodels.stats.multitest import fdrcorrection


def perform_de(adata_fg, adata_bg, tests=None):
    if tests is None:
        tests = ['ranksums']


    x = adata_fg.X
    if issparse(x):
        x = x.toarray()
    y = adata_bg.X
    if issparse(y):
        y = y.toarray()

    mean_fg = np.mean(x, axis=0)
    mean_bg = np.mean(y, axis=0)

    fold_change = mean_fg / mean_bg
    log2fc = np.log2(fold_change)

    nz_fg = np.count_nonzero(x, axis=0)
    nz_bg = np.count_nonzero(y, axis=0)

    pct_expr_fg = (nz_fg / adata_fg.X.shape[0]) * 100
    pct_expr_bg = (nz_bg / adata_bg.X.shape[0]) * 100

    delta_pct_expr = pct_expr_fg - pct_expr_bg

    de = pd.DataFrame(np.reshape(log2fc, (len(adata_fg.var_names), 1)), columns=['log2fc'], index=list(adata_fg.var_names))
    for test in tests:
        if test == 'ranksums':
            rv = sps.ranksums(x, y, axis=0)
            p = rv.pvalue
            _, pfdr = fdrcorrection(p)
            de[f'{test}-p'] = p
            de[f'{test}-fdr-p'] = pfdr
        elif test == 'ttest':
            rv = sps.ttest_ind(x, y, axis=0)
            p = rv.pvalue
            _, pfdr = fdrcorrection(p)
            de[f'{test}-p'] = p
            de[f'{test}-fdr-p'] = pfdr
        else:
            raise Exception(f'statistical test {test} not supported')

    de['fg_lin-avg'] = mean_fg
    de['bg_lin-avg'] = mean_bg
    de['fg_pct-expr'] = pct_expr_fg
    de['bg_pct-expr'] = pct_expr_bg
    de['delta-pct-expr'] = delta_pct_expr
    return de


def differential_expression(adata, comparisons, **de_kwargs):
    dedict = {}
    for cname, (fgf, bgf) in comparisons.items():
        adata_fg = adata_filter(adata, fgf)
        adata_bg = adata_filter(adata, bgf)
        dedict[cname] = perform_de(adata_fg, adata_bg, **de_kwargs)
    return dedict


def flag_de(de_tables, log2fc_thresh=1.0, p_column='ranksums-fdr-p', p_thresh=.05):
    if not isinstance(de_tables, dict):
        de_tables = {'single': de_tables}

    def is_de(row):
        if np.abs(row['log2fc']) > log2fc_thresh and row[p_column] < p_thresh:
            return True
        return False

    for de in de_tables.values():
        de['is-de'] = [is_de(row) for i, row in de.iterrows()]


def get_de_genes(de_tables, genesets=None, genes=None, only_up=False):
    """
    """
    # usually takes a dict of data frames, but could take a single just as well
    if isinstance(de_tables, pd.DataFrame):
        de_tables = {'single': de_tables}

    if genes is None:
        genes = []

    if genesets:
        for gsl, gsn in genesets.items():
            gs = get_geneset(gsn)
            genes.extend(gs)
    genes = set(genes)

    de_genes = set()

    for de_name, de in sorted(de_tables.items(), key=lambda x: x[0]):
        if 'is-de' not in de:
            raise Exception(f'must use flag_de function to add is-de column to data')
        de = de[de['is-de']]
        if only_up:
            de = de[de['log2fc'] > 0]
        for g in de.index:
            if len(genes) > 0:
                if g in genes:
                    de_genes.add(g)
            else:
                de_genes.add(g)
    return de_genes

def summarize_de_genes(de_tables, genesets=None, genes=None, always_include_manual_genes=True, only_up=True):
    by_geneset = {}
    all_de_genes = set()
    for gsl, gsn in genesets.items():
        de_genes = get_de_genes(de_tables, genesets={gsl: gsn}, only_up=only_up)
        by_geneset[gsl] = de_genes
        all_de_genes = all_de_genes.union(de_genes)

    if genes:
        de_genes = get_de_genes(de_tables, genes=genes, only_up=only_up)
        by_geneset['manual'] = de_genes
        all_de_genes = all_de_genes.union(de_genes)
        if always_include_manual_genes:
            by_geneset['manual'] = set(genes)
            all_de_genes = all_de_genes.union(set(genes))

    rows = []
    columns_start = ['index', 'gene', 'gene_group', 'cell_group', 'log2fc']
    columns_mid = None
    columns_end = ['is-de']
    for de_name, de in de_tables.items():
        cell_group = de_name
        if ':' in cell_group:
            cell_group = cell_group[:cell_group.index(':')]
        for g, row in de.iterrows():
            if g not in all_de_genes:
                continue
            geneset = None
            for gsl, gs in by_geneset.items():
                if g in gs:
                    geneset = gsl
                    break

            out_row = [g, g, geneset, cell_group, row['log2fc']]
            pheaders = []
            for c in [c for c in de.columns if c.endswith('-p')]:
                pheaders.append(c)
                out_row.append(row[c])
            if columns_mid is None:
                columns_mid = pheaders
            out_row.append(row['is-de'])
            rows.append(out_row)

    df = pd.DataFrame(rows, columns=columns_start + columns_mid + columns_end)
    df.set_index('index', inplace=True)
    df.sort_values(by=['gene_group', 'gene', 'cell_group'], inplace=True)
    return df


