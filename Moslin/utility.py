import warnings
from typing import List, Literal, Optional, Tuple
import warnings

import moscot as mt
import moscot.plotting as mtp
import ot
import tqdm
from moscot.problems.time import LineageProblem, TemporalProblem

import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import euclidean_distances

import matplotlib.pyplot as plt

import scanpy as sc
from anndata import AnnData

def compute_w2(
    p_1: np.ndarray,
    p_2: np.ndarray,
    cost: np.ndarray,
    metric_type: Literal["descendant", "ancestor"],
    scale_by_marginal: bool = True,
) -> List[float]:
    
    '''Iterate over the rows or columns in two coupling matrices and compute the earth mover’s distance (EMD).'''
    errors = []

    # transpose for ancestors
    if metric_type == "ancestor":
        p_1 = p_1.T
        p_2 = p_2.T

    # iterate over all cells
    for i in tqdm.tqdm(range(p_1.shape[0])):
        # normalize to get distributions
        marginal = p_1[i].sum()
        dist_1 = p_1[i] / marginal
        dist_2 = p_2[i] / p_2[i].sum()

        # compute the EMD distance between the ground-truth and the prediction
        error, log = ot.emd2(
            dist_1.astype(float), dist_2.astype(float), cost.astype(float), log=True
        )

        # append the EMD (potentially scaled by the marignal)
        if log["warning"] is None:
            if scale_by_marginal:
                errors.append(marginal * error)
            else:
                errors.append(error)
        else:
            errors.append(np.nan)

    return errors

def get_cost(
    adata: AnnData,
    early_tp: float = 8,
    late_tp: float = 12,
    time_key: str = "time",
    rep: str = "X_pca",
    n_dim: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    
    '''Obtain RNA-based cost matrices at early and late time points from an AnnData object.'''
    # retrieve cost matrices
    embedding = np.asarray(adata.obsm[rep], dtype=float)

    # make sure we have enough dimensions
    if n_dim is None:
        n_dim = embedding.shape[1]
    elif embedding.shape[1] < n_dim:
        raise ValueError(f"You only have {embedding.shape[1]} dimensions in {rep}.")
    embedding = embedding[:, :n_dim]

    # restrict to early and late cells
    mask_early, mask_late = (
        adata.obs[time_key] == early_tp,
        adata.obs[time_key] == late_tp,
    )
    embedding_early, embedding_late = embedding[mask_early], embedding[mask_late]

    # compute the corresponding distance matrices
    cost_early = euclidean_distances(embedding_early, squared=True)
    cost_late = euclidean_distances(embedding_late, squared=True)

    return cost_early, cost_late

def compute_errors(
    adata: AnnData,
    p_1: np.ndarray,
    p_2: np.ndarray,
    early_tp: float = 8,
    late_tp: float = 12,
    time_key: str = "time",
    rep: str = "X_pca",
    n_dim: Optional[int] = None,
    scale_by_marginal: bool = True,
) -> Tuple[List[float], List[float]]:
    
    '''Compute ancestor and descendant error between two couplings, where one typically corresponds to the ground-truth.'''
    # check that the shapes match
    assert p_1.shape == p_2.shape

    # retrieve early and late cost matrices
    cost_early, cost_late = get_cost(
        adata=adata,
        early_tp=early_tp,
        late_tp=late_tp,
        time_key=time_key,
        rep=rep,
        n_dim=n_dim,
    )

    # compute the ancestor errors
    ancestor_errors = compute_w2(
        p_1=p_1,
        p_2=p_2,
        cost=cost_early,
        metric_type="ancestor",
        scale_by_marginal=scale_by_marginal,
    )

    # compute the descendant errors
    descendant_errors = compute_w2(
        p_1=p_1,
        p_2=p_2,
        cost=cost_late,
        metric_type="descendant",
        scale_by_marginal=scale_by_marginal,
    )

    return ancestor_errors, descendant_errors

def compute_pearsonNorm_pca_umap_leiden(
        adata, 
        n_highly_variable,
        n_pc, 
):
    adata.layers['raw'] = adata.X.copy()
    sc.experimental.pp.recipe_pearson_residuals(adata, n_comps=50, n_top_genes=2000, chunksize=1000, clip=np.sqrt(adata.n_obs))


def matrixfy_character_obsm(
        adata,
        key: str | None,
        lineage_barcode_obsm_name: str | None = 'characters',
        special_symbols_characterization: dict | None = {"*": 0, "-": -1, "!": 999}, 
):
    '''
    Processing lineage_barcode to remove any non numerical values in the matrix. 
    Note: Special symbols meaning was estimated. If actually annotation exist (from pipeline documentation, for example),
    should pass in that instead. 

    Parameters
    ----------
    adata 
        An adata containing gene and lineage information from ALL days, and especially should contain an obsm for lineage barcode character matrix.

    key
        A name for the modified obsm matrix to be stored in. By default replaces the current character matrix. 

    lienage_barcode_obsm_name
        Name of the character matrix obsm
    
    special_symbols_characterization
        A dictionary containing all the symbols with their respecive numerical value to replace. By default, replaces the 
        *, -, ! to 0, -1, and 999 respectively. 
        Meaning of symbols:
        "*": 0, no change
        "-": -1, not detected
        "!": 999, detected, but not editable
    
    Returns
    -------
    adata
        Adata with the updated character matrix, where it should only be compose of int64 datatype; stored in `key` if specified, 
        otherwise replace the existing character matrix 
    '''
    if key == None:
        key = lineage_barcode_obsm_name

    # Extracting the character matrix 
    chars = adata.obsm[lineage_barcode_obsm_name]

    for key, value in special_symbols_characterization.items():
        chars = np.where(chars == key, value, chars)
    
    adata.obsm[key] = chars.astype(int)

    if (adata.obsm[key].dtype) != np.int32:
        warnings.warn("Matrix is not fully transform to numerical values. There might be a symbol that was not replaced.")

    return adata
