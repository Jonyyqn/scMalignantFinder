import scanpy as sc
from pyscenic.aucell import aucell, derive_auc_threshold
from ctxcore.genesig import GeneSignature

def aucell_cal(adata, gmt, norm_type=False, adapt_signatures=True):
    """
    Perform AUCell-based enrichment scoring on an AnnData object.

    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix.
    gmt : str
        Path to gene set in GMT format.
    norm_type : bool, optional
        Whether to normalize total counts (default: False).
    adapt_signatures : bool, optional
        Whether to filter gene sets to those with sufficient overlap (default: True).

    Returns:
    --------
    adata : AnnData
        AnnData object with AUCell scores added to `.obs`.
    """
    if norm_type:
        sc.pp.normalize_total(adata, target_sum=1e4)

    if adapt_signatures:
        geneset_dict = {}
        with open(gmt) as f:
            for line in f:
                data = line.strip().split('\t')
                geneset_dict[data[0]] = data[2:]

        gmt_tmp = f'{gmt}.tmp'
        with open(gmt_tmp, 'w') as f:
            for title, genes in geneset_dict.items():
                overlap_genes = list(set(genes) & set(adata.var_names))
                if len(overlap_genes) / len(genes) >= 0.5:
                    f.write('{0}\t{0}\t{1}\n'.format(title, '\t'.join(overlap_genes)))

        gs = GeneSignature.from_gmt(gmt_tmp, field_separator="\t", gene_separator="\t")
    else:
        gs = GeneSignature.from_gmt(gmt, field_separator="\t", gene_separator="\t")

    df = adata.to_df()
    percentiles = derive_auc_threshold(df)

    scores = aucell(
        exp_mtx=df,
        signatures=gs,
        auc_threshold=percentiles[0.01],
        seed=2,
        normalize=True,
        num_workers=1
    )

    for col in scores.columns:
        adata.obs[col] = scores.loc[adata.obs_names, col].tolist()

    return adata