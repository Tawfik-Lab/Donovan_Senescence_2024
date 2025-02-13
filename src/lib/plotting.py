import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import issparse
import plotly.graph_objects as go
import plotly.express as px
from util import get_geneset
from collections import OrderedDict


def plot_de_genes(adata, obs, de_tables, genesets=None, genes=None, sort_genes_by='change', **kwargs):
    """
    """

    var_names = set(adata.var_names)
    gss = OrderedDict()
    gsa = set()
    if genesets:
        for gsl, gsn in genesets.items():
            try:
                gs = get_geneset(gsn)
            except:
                raise(f'Unable to find geneset {gsn}')
            gss[gsl] = gs
            gsa = gsa.union(set(gs))


    def superset(dos):
        ss = set()
        for k, vs in dos.items():
            for v in vs:
                ss.add(v)
        return ss

    gk = OrderedDict()
    max_expr = dict()

    data_rows = []
    data_headers_start = ['index', 'gene', 'gene_group', 'cell_group', 'log2fc']
    data_headers_mid = None
    data_headers_end = ['is-de']

    for de_name, de in sorted(de_tables.items(), key=lambda x: x[0]):
        de = de[de['is-de']]
        de = de[de['log2fc'] > 0.0]
        de = de.sort_values(by='log2fc', ascending=False)

        for g, row in de.iterrows():

            for gsl, gs in gss.items():
                if g in gs and ( (not genes) or (g not in genes) ) and (g not in superset(gk)):
                    if gsl not in gk:
                        gk[gsl] = set()
                    gk[gsl].add(g)
                    if sort_genes_by == 'change':
                        gv = row['log2fc']
                    elif sort_genes_by == 'expr':
                        gv = row['fg_lin-avg']
                    else:
                        continue

                    if g not in max_expr or gv > max_expr[g]:
                        max_expr[g] = gv

    plot_genes = []
    var_group_positions = []
    var_group_labels = []

    if genes is not None:
        plot_genes = [x for x in genes]

    i = len(plot_genes)
    for gsl, gs in gk.items():
        if len(gs) < 1:
            continue
        if sort_genes_by in ['change', 'expr']:
            gs = sorted(gs, key=lambda g: max_expr[g], reverse=True)
        else:
            gs = sorted(gs)

        vgp = [i]
        for g in gs:
            i += 1
            plot_genes.append(g)
        vgp.append(i - 1)
        var_group_positions.append(vgp)
        var_group_labels.append(gsl)

    return dotplot(
        adata,
        plot_genes,
        groupby=obs,
        max_scale='var',
        var_group_positions=var_group_positions,
        var_group_labels=var_group_labels,
        **kwargs,
    )

def dotplot(adata, var_names, max_scale=None, *args, **kwargs):
    groupby = kwargs.get('groupby', None)
    if groupby:
        if isinstance(groupby, str):
            groupby = [groupby]
        elif isinstance(groupby, list):
            if not all(isinstance(g, str) for g in groupby):
                raise ValueError("all elements in groupby must be strings")
        else:
            raise ValueError("groupby must be a string or a list of strings")
        
        for g in groupby:
            if g in adata.obs:
                adata.obs[g] = adata.obs[g].astype('category')
            else:
                raise ValueError(f"'{g}' is not a valid column in adata.obs")
                
    if max_scale == 'var':
        adata = adata[:, [x in set(var_names) for x in adata.var_names]].copy()
        X = adata.X
        if issparse(X):
            X = X.toarray()
        df = pd.DataFrame(X, columns=adata.var_names, index=adata.obs_names)
        df['group'] = adata.obs[kwargs['groupby']]
        gmin = df.groupby('group').mean().min(axis=0)
        div = df.groupby('group').mean().max(axis=0)
        adata.X = (X - gmin.values) / (div.values - gmin.values)
    elif max_scale == 'obs':
        adata = adata[:, [x in set(var_names) for x in adata.var_names]].copy()
        X = adata.X
        if issparse(X):
            X = X.toarray()
        df = pd.DataFrame(X, columns=adata.var_names, index=adata.obs_names)
        df['group'] = adata.obs[kwargs['groupby']]
        dd = df.groupby('group').mean().max(axis=1).to_dict()
        gmin = df.groupby('group').mean().min(axis=0).to_dict()
        div = [dd[g] for g in adata.obs[kwargs['groupby']]]
        div = np.reshape(np.array(div), (adata.shape[0], 1))
        gmin = [gmin[g] for g in adata.obs[kwargs['groupby']]]
        gmin = np.reshape(np.array(gmin), (adata.shape[0], 1))
        adata.X = (X - gmin) / (div - gmin)
    elif max_scale is not None:
        raise Exception(f'invalid value {max_scale} for max_scale: should be obs or var')
        
    return sc.pl.dotplot(adata, var_names, *args, **kwargs)


def plot_score_heatmap(adata, xaxis_key, yaxis_key, score_key, y_as=None, x_as=None, layout_kwargs=None):
    rows = []
    xs = list(sorted(adata.obs[xaxis_key].unique()))
    ys = list(sorted(adata.obs[yaxis_key].unique()))

    if x_as == 'int':
        xs = list(sorted(xs, key=lambda z: int(z)))
    if y_as == 'int':
        ys = list(sorted(ys, key=lambda z: int(z)))
        
    for yk in ys:
        row = []
        for xk in xs:
            adata_xy = adata[adata.obs[xaxis_key] == xk].copy()
            adata_xy = adata_xy[adata_xy.obs[yaxis_key] == yk]
            value = adata_xy.obs[score_key].mean()
            row.append(value)
        rows.append(row)
    heatmap_values = pd.DataFrame(np.round(np.array(rows), 2), index=ys, columns=xs)

    fig = px.imshow(
        heatmap_values,
        text_auto=True, 
        color_continuous_scale='RdBu_r',
        color_continuous_midpoint=0.0,
    )
    title = f'Mean {score_key} Score'
    fig.update_layout(title=title, height=800, width=800, font_size=18)
    fig.update_xaxes(title=xaxis_key, tickangle=-90)
    fig.update_yaxes(title=yaxis_key)
    if layout_kwargs:
        fig.update_layout(**layout_kwargs)
    return fig

