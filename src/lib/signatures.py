from scipy.sparse import issparse
import pandas as pd
import numpy as np


def parse_signature(sig):
    sp = sig.split(',')
    parsed = []
    for element in sp:
        if element.startswith('-'):
            parsed.append(('-', element[1:]))
        else:
            parsed.append(('+', element))
    return parsed

def assess_signatures(adata, sigs):
    for sig in sigs:
        parsed = parse_signature(sig)
        genes = {p[1]: p[0] for p in parsed}
        adata_sig = adata[:, [x in genes for x in adata.var_names]].copy()
        if adata_sig.shape[1] != len(parsed):
            raise Exception(f'not all genes found, only have {list(adata_sig.var_names)} out of {list(genes.keys())}')
        if issparse(adata_sig.X):
            adata_sig.X = adata_sig.X.toarray()
        sigdata = pd.DataFrame(adata_sig.X, index=list(adata_sig.obs_names), columns=list(adata_sig.var_names))
        sign_order = []


        # compute additive
        sigdata_ad = sigdata.copy()
        for c in list(sigdata_ad.columns):
            sign = genes[c]
            if sign == '-':
                sigdata_ad[c] = -1 * sigdata_ad[c]
                sign_order.append(0)
            else:
                sign_order.append(1)
        print(f'{sig} {parsed} {list(sigdata.columns)} {sign_order}')

        additive = sigdata_ad.sum(axis=1)

        # compute binary 
        sigdata_binary = (sigdata > 0.0).astype('int')
        sigdata_binary_xor = np.logical_xor(sigdata_binary, sign_order)
        binary = 1.0 - (sigdata_binary_xor.sum(axis=1) / sigdata_binary.shape[1])

        # compute exclusive
        exclusive = (binary == 1.0).astype('int') * additive

        adata.obs[f'ex:{sig}'] = exclusive
        adata.obs[f'ad:{sig}'] = additive
        adata.obs[f'bi:{sig}'] = binary

def summarize_signatures(adata_or_obs, gene_signatures, groupby, include_columns=None, totals_groupby=None, sort_by=None, sort_ascending=None):
    if include_columns is None:
        include_columns = []

    if not isinstance(groupby, list):
        groupby = [groupby]

    if isinstance(adata_or_obs, pd.DataFrame):
        obs = adata_or_obs
    else:
        obs = adata_or_obs.obs

    sig_columns = [f'ex:{sig}' for sig in gene_signatures]

    sig_df = obs[groupby + include_columns + sig_columns].copy()
    
    # aggregations configuration
    agg = {}
    agg['total'] = pd.NamedAgg(column=groupby[0], aggfunc='count')
    for sig in gene_signatures:
        agg[f'{sig}_num_cells_positive'] = pd.NamedAgg(column=f'ex:{sig}', aggfunc=lambda x: (x > 0.0).astype('int').sum())
    for col in include_columns:
        dt = sig_df[col].dtype
        if dt == 'float' or dt == 'int':
            agg[col] = pd.NamedAgg(column=col, aggfunc=lambda x: x.mean())
        else:
            agg[col] = pd.NamedAgg(column=col, aggfunc=lambda x: ','.join(sorted(set(x))))

    sig_df = sig_df.groupby(groupby).agg(**agg)

    # get the indexed columns back into the data
    sig_df_no_i = sig_df.reset_index()
    sig_df_no_i.index = sig_df.index
    sig_df = sig_df_no_i


    # if we want to calculate totals from some subset of groupby
    totals_per_total_group = None
    total = None
    if totals_groupby is not None:
        if totals_groupby not in groupby:
            raise Exception(f'totals_groupby must be one of the initial groupby: {totals_groupby} not in {groupby}')
        totals_per_total_group = sig_df[['total']]
        totals_per_total_group.reset_index(inplace=True)
        totals_per_total_group = totals_per_total_group.groupby(totals_groupby).agg({'total': 'sum'}).to_dict()['total']
    else:
        total = sig_df['total'].sum()

    for sig in gene_signatures:
        sig_df[f'{sig}_percent_group_positive'] = (sig_df[f'{sig}_num_cells_positive'] / sig_df['total']) * 100
        if totals_groupby and totals_per_total_group:
            sig_df[f'{sig}_percent_all_positive'] = [(p / totals_per_total_group[tg]) * 100 for (p, tg) in zip(sig_df[f'{sig}_num_cells_positive'], sig_df[totals_groupby])]
        else:
            sig_df[f'{sig}_percent_all_positive'] = [(p / total) * 100 for p in sig_df[f'{sig}_num_cells_positive']]

    first_columns = groupby + include_columns + ['total']
    sig_columns = []
    for sig in gene_signatures:
        sig_columns.extend(list(sorted([c for c in sig_df.columns if c.startswith(f'{sig}_')])))

    if sort_by:
        if not isinstance(sort_by, list):
            sort_by = [sort_by]
        sort_by_call = []
        for s in sort_by:
            sn = s
            if s in groupby:
                sn = f'{s}_sort'
                sig_df[sn] = sig_df[s]
            sort_by_call.append(sn)
        sort_kwargs = {}
        if sort_ascending:
            sort_kwargs['ascending'] = sort_ascending
        sig_df.sort_values(by=sort_by_call, inplace=True, **sort_kwargs)


    sig_df = sig_df[first_columns + sig_columns]
    return sig_df
