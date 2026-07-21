import gc
import numpy as np
import pandas as pd
import networkx as nx
from pathlib import Path
import pycea as py
import sys
sys.path.append("/project/imoskowitz/yubin/Lineage_Tree_Construction")
from lineage_utilities.tree_functions import *

def analyze_tree(clone_tdata, clone_key, color, distance_metric='lca', n_permutations=500, random_state=0, ):
    """
    Single-pass per-tree analysis: runs shared preprocessing once, then
    computes leaf-depth records, polytomy/resolution stats, and a
    same-type-vs-different-type LCA permutation test — all reusing the
    same tree_distance/LCA matrix and depth annotation.

    Parameters
    ----------
    color 
        What to color the leaves as well as conduct comparison by. Usually the celltype or some other sort of cluster labeling. 
    Returns
    -------
    dict with keys: "leaf_depth_rows", "resolution_stats", "permutation_stats"
    """
    # --- shared preprocessing (once per tree) ---
    add_obs_to_tree(clone_tdata, keys=[color, "germ_layer"]) # Not sure what germ_layer is for, maybe remove?
    py.pp.add_depth(clone_tdata)
    py.tl.tree_distance(clone_tdata, tree=clone_key, metric=distance_metric, key_added=distance_metric)

    tree_graph = clone_tdata.obst[clone_key]
    leaf_order = py.get.leaves(clone_tdata, tree=clone_key)
    n = len(leaf_order)

    leaf_idx_map = {name: i for i, name in enumerate(clone_tdata.obs_names)}
    leaf_indices = np.array([leaf_idx_map[l] for l in leaf_order])

    cell_types = clone_tdata.obs.loc[leaf_order, color].astype(str).to_numpy()
    valid_mask_leaf = ~np.isin(cell_types, ["nan", "None", "", "<NA>"])

    # --- leaf depth table (question 1) ---
    depths = clone_tdata.obs.loc[leaf_order, "depth"].to_numpy()
    max_depth = depths.max() if len(depths) else 1.0
    leaf_depth_rows = [
        {"tree": clone_key, "cell_type": ct, "depth": d,
         "depth_norm": d / max_depth if max_depth > 0 else 0.0, "n_leaves": n}
        for leaf, d, ct, valid in zip(leaf_order, depths, cell_types, valid_mask_leaf) if valid
    ]

    # --- polytomy / resolution diagnostic (question 3) ---
    out_degrees = dict(tree_graph.out_degree())
    polytomy_nodes = {node: d for node, d in out_degrees.items() if d > 2}
    node_depths = nx.get_node_attributes(tree_graph, "depth")
    polytomy_depths = [node_depths[node] for node in polytomy_nodes if node in node_depths]
    n_internal = sum(1 for node, d in out_degrees.items() if d > 0)

    resolution_stats = {
        "tree": clone_key,
        "n_leaves": n,
        "n_polytomies": len(polytomy_nodes),
        "n_internal_nodes": n_internal,
        "polytomy_fraction": len(polytomy_nodes) / max(1, n_internal),
        "mean_polytomy_depth_norm": (
            float(np.mean(polytomy_depths) / max_depth) if polytomy_depths and max_depth > 0 else np.nan
        ),
    }

    # --- pairwise LCA + permutation test (question 2) ---
    lca_depths = clone_tdata.obsp["lca_distances"]
    is_sparse = hasattr(lca_depths, "toarray")
    if is_sparse:
        lca_depths = lca_depths.tocsr()

    sub = lca_depths[leaf_indices, :][:, leaf_indices]
    if is_sparse:
        sub = sub.toarray()

    tri_i, tri_j = np.triu_indices(n, k=1)
    pair_depths = sub[tri_i, tri_j]
    valid_pair = valid_mask_leaf[tri_i] & valid_mask_leaf[tri_j]

    same_type = cell_types[tri_i] == cell_types[tri_j]
    same_mean = pair_depths[same_type & valid_pair].mean() if (same_type & valid_pair).any() else np.nan
    diff_mean = pair_depths[~same_type & valid_pair].mean() if (~same_type & valid_pair).any() else np.nan
    observed_diff = same_mean - diff_mean

    rng = np.random.default_rng(random_state)
    null_diffs = np.empty(n_permutations)
    for p in range(n_permutations):
        shuffled = rng.permutation(cell_types)
        same_type_p = shuffled[tri_i] == shuffled[tri_j]
        valid_p = ~np.isin(shuffled[tri_i], ["nan", "None", "", "<NA>"]) & \
                  ~np.isin(shuffled[tri_j], ["nan", "None", "", "<NA>"]) 
        same_mean_p = pair_depths[same_type_p & valid_p].mean() if (same_type_p & valid_p).any() else np.nan
        diff_mean_p = pair_depths[~same_type_p & valid_p].mean() if (~same_type_p & valid_p).any() else np.nan
        null_diffs[p] = same_mean_p - diff_mean_p

    p_value = (np.nansum(null_diffs >= observed_diff) + 1) / (n_permutations + 1)

    permutation_stats = {
        "tree": clone_key,
        "same_type_mean_lca_depth": same_mean,
        "diff_type_mean_lca_depth": diff_mean,
        "observed_diff": observed_diff,
        "p_value": p_value,
    }

    # free the dense matrix explicitly before returning
    del sub, lca_depths

    return {
        "leaf_depth_rows": leaf_depth_rows,
        "resolution_stats": resolution_stats,
        "permutation_stats": permutation_stats,
    }


def run_full_tree_analysis(adata, tree_key='tree', celltype_key='cardiac_labels',
                             output_dir=None, n_permutations=500, excluded_trees=None):
    """
    Loops over every tree once, running all diagnostics, and writes three
    combined output tables incrementally. Trees in `excluded_trees` are
    skipped entirely (e.g. very large trees to handle separately later).

    Parameter
    ---------
    adata 
        tdata that contains obst for clones
    
    tree_key
        the obs column that contain the name for tree
    
    
    """
    print("starting function")
    excluded_trees = excluded_trees or set()
    leaf_depth_csv = output_dir + "leaf_depth_all_trees.csv"
    resolution_csv = output_dir + "resolution_all_trees.csv"
    permutation_csv = output_dir + "permutation_all_trees.csv"

    for path in [leaf_depth_csv, resolution_csv, permutation_csv]:
        p = Path(path)
        if p.exists():
            p.unlink()

    tree_sizes = adata.obs[tree_key].value_counts()
    failed_trees = []
    skipped_trees = []

    for clone_key in adata.obs[tree_key].unique():
        if pd.isna(clone_key) or clone_key == 'nan':
            continue

        if clone_key in excluded_trees:
            print(f"skipping {clone_key} (excluded)")
            skipped_trees.append(clone_key)
            continue

        n_cells = tree_sizes.get(clone_key, 0)
        print(f"analyzing {clone_key} ({n_cells} cells)")

        try:
            clone_tdata = adata[adata.obs[tree_key] == clone_key].copy()
            result = analyze_tree(clone_tdata, clone_key, color=celltype_key,
                                   n_permutations=n_permutations)

            pd.DataFrame(result["leaf_depth_rows"]).to_csv(
                leaf_depth_csv, mode='a', header=not Path(leaf_depth_csv).exists(), index=False)
            pd.DataFrame([result["resolution_stats"]]).to_csv(
                resolution_csv, mode='a', header=not Path(resolution_csv).exists(), index=False)
            pd.DataFrame([result["permutation_stats"]]).to_csv(
                permutation_csv, mode='a', header=not Path(permutation_csv).exists(), index=False)

            del clone_tdata, result

        except Exception as e:
            print(f"ERROR on {clone_key}, skipping: {e}")
            failed_trees.append(clone_key)

        finally:
            gc.collect()

        print(f"complete {clone_key}")

    print("Skipped trees (excluded):", skipped_trees)
    print("Failed trees:", failed_trees)
    return {
        "leaf_depth_csv": leaf_depth_csv,
        "resolution_csv": resolution_csv,
        "permutation_csv": permutation_csv,
        "failed_trees": failed_trees,
        "skipped_trees": skipped_trees,
    }

adata = td.read_h5td("/project/imoskowitz/yubin/Lineage_Tree_Construction/output_data/Processed_data/E8_5.h5td")
data_dir = "output_data"
plot_dir = "output_plot"
base_path = "/project/imoskowitz/yubin/Lineage_Tree_Construction/"
output_path_data = base_path+data_dir+"/Trees/Robin_Pijuan_celltypist/Manual_Cardiac_Annotation/E8_5/"
output_path_plot = base_path+plot_dir+"/Trees/Robin_Pijuan_celltypist/Manual_Cardiac_Annotation/E8_5/"


# EXCLUDED_TREES = {"E8.5-R3-C1", "E8.5-R2-C1", "E8.5-R2-C2"} Run one excluded these trees
# Run two will compute these trees but exclude the rest that are already computed
EXCLUDED_TREES = set(adata.obs['tree']) - {"E8.5-R3-C1", "E8.5-R2-C1", "E8.5-R2-C2"} 
paths = run_full_tree_analysis(
    adata,
    tree_key='tree',
    celltype_key='cardiac_labels',
    output_dir=output_path_plot,
    n_permutations=500,
    excluded_trees=EXCLUDED_TREES
)