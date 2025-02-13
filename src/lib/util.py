import os
import os.path
import json


def adata_filter_mask(adata, filter_dict):
    """
    Filter dict is a dictionary of obs keys to 'filter values'.
    Filter values can be 
    * a scalar, in which case the obs value must equal it.
    * a list in which case the obs value must be in that list
    * a callable, which will be called with the obs value and returns True or False
    """

    mask = []
    rn = 0
    vicache = {}
    for i, row in adata.obs.iterrows():
        rv = True
        for fk, fv in filter_dict.items():
            if fk in row:
                row_value = row[fk]
            elif fk in adata.var_names:
                if fk in vicache:
                    vi = vicache[fk]
                else:
                    vi = list(adata.var_names).index(fk)
                    vicache[fk] = vi
                row_value = adata.X[rn,vi]
            else:
                raise Exception(f'filter key {fk} not in adata')

            if callable(fv):
                if not fv(row_value):
                    rv = False
                    break
            elif isinstance(fv, list):
                if row_value not in fv:
                    rv = False
                    break
            elif isinstance(fv, dict):
                for fvk, fvv in fv.items():
                    if fvk == 'gt':
                        if not row_value > fvv:
                            rv = False
                    # could contain other options
                    else:
                        raise Exception(f'unsupported filter component {fv} in filter {filter_dict}')
            elif row_value != fv:
                rv = False
                break
        mask.append(rv)
        rn += 1
    return mask

def adata_filter(adata, filter_dict):
    """
    Filter dict is a dictionary of obs keys to 'filter values'.
    Filter values can be 
    * a scalar, in which case the obs value must equal it.
    * a list in which case the obs value must be in that list
    * a callable, which will be called with the obs value and returns True or False
    """

    mask = adata_filter_mask(adata, filter_dict)
    return adata[mask]


def get_geneset(geneset_name):
    dirname = os.path.dirname(os.path.abspath(__file__))
    geneset_dirname = os.path.abspath(os.path.join(dirname, '../../input_data'))
    geneset_filename = os.path.join(geneset_dirname, f'{geneset_name}.json')
    with open(geneset_filename, 'r') as fh:
        geneset = json.load(fh)
    return geneset
