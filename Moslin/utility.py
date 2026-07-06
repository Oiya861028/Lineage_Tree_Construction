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
from fastrna.core import fastrna_hvg, fastrna_pca

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

def compute_hvg_pca_fastRNA(
        adata, 
        batch_key,
        n_highly_variable = 2000,
):
    """
    Computes the hvg and pca using fastRNA, a more memory and speed efficient method than conventional methods. Takes in raw data.

    Parameters
    ----------
    adata
        Contains raw data in the X. 

    batch_key
        Name of obs that labels the cells' batch
    
    n_highly_variable
        Number of highly variable genes to use for pca

    return
    ------
    adata 
        Containing X_pca and pc variance.
    """

    gene_exp = adata.X

    # Getting index of cells by batch
    batch_label = pd.factorize(adata.obs[batch_key])[0] # change the column name (Method) if required,
    batch_label
    idx_sort = batch_label.argsort()

    # Sorting cells by batch
    batch_label_sort = batch_label[idx_sort]
    gene_exp_sort = gene_exp[idx_sort,:].T # note: genes x cells, so sort columns

    # Compute hvg
    gene_vars = fastrna_hvg(gene_exp_sort, batch_label_sort)
    gene_idx_var = gene_vars.argsort()[::-1]
    gene_exp_hvg = gene_exp_sort[gene_idx_var[:n_highly_variable], :]  
    gene_exp_hvg.sort_indices() # required after gene selection

    # Compute pca
    numi = np.asarray(gene_exp_sort.sum(axis=0)).ravel() # size factors (total UMI per cell)
    eig_val, eig_vec, pca, rrt = fastrna_pca(gene_exp_hvg, numi, batch_label_sort)

    # Storing Results

    # Storing Pca
    idx_unsort = np.argsort(idx_sort) # Need to undo sorting from earlier
    adata
    adata.obsm["X_pca"] = pca[idx_unsort, :]

    # Store variance explained (eigenvalues) — scanpy expects ratio and ratio_cumsum
    total_variance = eig_val.sum()
    adata.uns["pca"] = {
        "variance": eig_val,
        "variance_ratio": eig_val / total_variance,
    }

    # Store gene loadings (eigenvectors) — shape should be (n_genes, n_pcs)
    # eig_vec columns correspond to PCs, rows correspond to the selected HVGs
    # Need to map back to full gene space
    loadings_full = np.zeros((adata.n_vars, eig_vec.shape[1]))
    loadings_full[gene_idx_var[:n_highly_variable], :] = eig_vec
    adata.varm["PCs"] = loadings_full

    return adata



def matrixfy_character_obsm(
        adata,
        key: str | None = None,
        lineage_barcode_obsm_name: str | None = 'X_characters',
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

    for k, v in special_symbols_characterization.items():
        chars = np.where(chars == k, v, chars)
    
    adata.obsm[key] = chars.astype(int)

    if (adata.obsm[key].dtype) != np.int32:
        warnings.warn("Matrix is not fully transform to numerical values. There might be a symbol that was not replaced.")

    return adata

def add_obs_to_tree(tdata, keys):
    """Manually push obs columns onto tree nodes, since 0.2.0 has no add_obs_annotation."""
    for tree_name, tree in tdata.obst.items():
        for node in tree.nodes:
            # leaf nodes are named by their cell barcode, matching tdata.obs.index
            if node in tdata.obs.index:
                for key in keys:
                    tree.nodes[node][key] = tdata.obs.loc[node, key]

def add_node_attrs_to_edges(tdata, keys):
    """Copy node attributes onto their outgoing edges, which is what pl.branches() reads."""
    for tree_name, tree in tdata.obst.items():
        for parent, child in tree.edges:
            for key in keys:
                # branches() colors by the child node's attribute
                val = tree.nodes[child].get(key)
                if val is not None:
                    tree.edges[parent, child][key] = val

def flag_mixed_nodes(tdata, tree_key, cell_type_key="cell_type"):
    """
    After ancestral_states, override cell_type to 'mixed' for any internal
    node whose direct children have more than one distinct cell type.
    Leaves are never overridden.
    """
    tree = tdata.obst[tree_key]
    for node in tree.nodes:
        if tree.out_degree(node) == 0:
            continue  # skip leaves
        child_types = set(
            tree.nodes[c].get(cell_type_key)
            for c in tree.successors(node)
            if tree.nodes[c].get(cell_type_key) is not None
        )
        if len(child_types) > 1:
            tree.nodes[node][cell_type_key] = "mixed"

def add_gene_expression_to_tree(tdata, tree_key, gene, layer=None):
    """
    For each leaf node, add mean expression of gene as a node attribute.
    For internal nodes, use the mean across all descendant leaves.

    Parameters
    ----------
    tdata : TreeData
    tree_key : str
    gene : str, must be in tdata.var_names
    layer : str or None, if None uses tdata.X
    """
    tree = tdata.obst[tree_key]

    # get expression vector for all cells
    gene_idx = tdata.var_names.get_loc(gene)
    if layer is not None:
        expr = np.array(tdata.layers[layer][:, gene_idx]).flatten()
    else:
        import scipy.sparse as sp
        X = tdata.X
        expr = np.array(X[:, gene_idx].todense()).flatten() if sp.issparse(X) else X[:, gene_idx]

    # map barcode -> expression value
    barcode_to_expr = dict(zip(tdata.obs.index, expr))

    # for each node, compute mean expression across all descendant leaves
    for node in reversed(list(nx.topological_sort(tree))):
        if tree.out_degree(node) == 0:
            # leaf: look up directly
            val = barcode_to_expr.get(node, np.nan)
        else:
            # internal: mean over all descendant leaves
            desc_leaves = [
                n for n in nx.descendants(tree, node)
                if tree.out_degree(n) == 0
            ]
            vals = [barcode_to_expr[l] for l in desc_leaves if l in barcode_to_expr]
            val = np.mean(vals) if vals else np.nan
        tree.nodes[node][gene] = val