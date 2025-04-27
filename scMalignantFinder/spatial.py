import numpy as np
import scanpy as sc
import squidpy as sq
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from scMalignantFinder import classifier
from scipy.cluster.hierarchy import linkage, cut_tree
from scipy.spatial.distance import pdist
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
        num_workers=10
    )

    for col in scores.columns:
        adata.obs[col] = scores.loc[adata.obs_names, col].tolist()

    return adata


def image_cal(adata, library_id=None, scale=None, img=None):
    """
    Extract image features and calculate normalized image-based scores.

    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix.
    library_id : str, optional
        Library ID of spatial data (default: first available).
    scale : float, optional
        Image scale factor (default: from `adata.uns`).
    img : sq.im.ImageContainer, optional
        Squidpy image container (default: loaded from adata).

    Returns:
    --------
    adata : AnnData
        AnnData object with normalized image feature score added to `.obs['image_score']`.
    """
    if not library_id:
        library_id = list(adata.uns['spatial'].keys())[0]

    if not scale:
        scale = adata.uns['spatial'][library_id]['scalefactors']['tissue_hires_scalef']

    if not img:
        img = sq.im.ImageContainer(
            adata.uns['spatial'][library_id]['images']['hires'],
            scale=scale,
            library_id=library_id
        )

    # Extract summary image features
    sq.im.calculate_image_features(
        adata,
        img,
        features="summary",
        key_added="image_summary",
        n_jobs=1,
        scale=scale
    )

    # Normalize image-based score
    scaler = MinMaxScaler()
    adata.obs['image_score'] = scaler.fit_transform(
        adata.obsm['image_summary']['summary_ch-0_quantile-0.5'].values.reshape(-1, 1)
    )[:, 0]

    return adata


def region_identification(
    adata,
    features=['malignancy_probability', 'Malignant_up', 'image_score'],
    nclus=3,
    define_feature='Malignant_up',
    spatial_nn=True
):
    """
    Identify tumor regions by clustering and optional spatial smoothing.

    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix.
    features : list of str
        List of features to use for clustering (must be ≥2).
    nclus : int
        Number of clusters to cut hierarchical tree into.
    define_feature : str
        The feature used to define the malignant cluster (highest average).
    spatial_nn : bool
        Whether to refine labels using spatial neighbors (default: True).

    Returns:
    --------
    adata : AnnData
        AnnData object with `.obs['region_prediction']` added ('Malignant' or 'Normal').
    """
    assert len(features) >= 2, "Feature list must contain at least 2 elements."

    # Hierarchical clustering
    matrix = adata.obs[features].values
    dist_matrix = pdist(matrix, metric='euclidean')
    linkage_matrix = linkage(dist_matrix, method='ward')
    clusters = pd.DataFrame({
        'cluster': cut_tree(linkage_matrix, n_clusters=nclus).flatten()
    }, index=adata.obs_names)
    adata.obs['cluster'] = pd.Categorical(clusters['cluster'])
    
    # Find tumor-like cluster
    # tumor_cluster = pd.concat([clusters, adata.obs], axis=1).groupby('cluster')[define_feature].mean().argmax()
    tumor_cluster = adata.obs.groupby('cluster')[define_feature].mean().argmax()
    adata.obs['region_prediction'] = pd.Categorical(
        clusters['cluster'].map(lambda x: 'Malignant' if x == tumor_cluster else 'Normal'),
        categories=['Normal', 'Malignant']
    )

    if spatial_nn:
        # Build spatial graph
        sq.gr.spatial_neighbors(adata, n_rings=1, coord_type="grid", n_neighs=6)

        raw_cluster = np.array(adata.obs['region_prediction'].values)
        spots = adata.obs_names.values
        new_cluster = []

        for cluster, spot in zip(raw_cluster, spots):
            nn_indices = np.where(np.array(
                adata.obsp["spatial_connectivities"][np.where(spots == spot)[0]].todense()
            )[0] == 1)[0]

            nn_cluster = raw_cluster[nn_indices]

            if len(nn_cluster) < 2:
                new_cluster.append(cluster)
            else:
                unique_values, counts = np.unique(nn_cluster, return_counts=True)
                max_count_index = np.argmax(counts)
                most_common_element = unique_values[max_count_index]

                if counts[max_count_index] / len(nn_cluster) > 0.5:
                    new_cluster.append(most_common_element)
                else:
                    new_cluster.append(cluster)

        adata.obs['region_prediction'] = pd.Categorical(new_cluster, categories=['Normal', 'Malignant'])

    return adata