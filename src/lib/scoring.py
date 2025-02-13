import scanpy as sc

def score_within_key(adata, geneset, obs_key, score_key='score'):
    score_lookup = {}
    for k in adata.obs[obs_key].unique():
        adata_k = adata[adata.obs[obs_key] == k].copy()
        sc.pp.log1p(adata_k)
        sc.pp.scale(adata_k, zero_center=True)
        sc.tl.score_genes(
            adata_k,
            geneset,
            ctrl_size=len(geneset),
            score_name=score_key,
        )
        for i, s in zip(adata_k.obs.index, adata_k.obs[score_key]):
            score_lookup[i] = s
    adata.obs[score_key] = [score_lookup[i] for i in adata.obs.index]
