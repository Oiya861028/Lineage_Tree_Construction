import warnings
from typing import List, Literal, Optional, Tuple
import warnings

import treedata as td
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import euclidean_distances

import matplotlib.pyplot as plt

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


def plot_tree_custom(clone_tdata,
                     clone_key,
                     color,
                     palette: None | str | list = None,
                     n_leaf=5,
                     output_path_plot: str | None = None
):
    """
    Given a tdata with obst, plots:
        - the tree labeling the leaf by specified color,
        - the branch with the color of the common descendent color and grayed when descendent color are not the same.
        - dotplot of the average depth of each color within n_leaf,
        - with a histogram summarizing the depth counts

    Parameters
    ----------
    clone_tdata (tdata)
        tdata containing lineage tree
    clone_key (str)
        name for the lineage tree
    color (str)
        obs name for the cell/leaf to be colored by
    palette (str | list | None)
        color palette for the obs. If None, attempt to find the default color in .uns
    n_leaf
        Number of leaf to cluster by to calculate the average depth of. This range will be iterated from the top to the bottom of the tree
    output_path_plot
        Saving plot on file location. If None, show plot instead.

    Returns
    -------
    dict with keys:
        "ax_tree"     : the tree axis (only meaningful if output_path_plot is None)
        "points"      : list of (depth, leaf_position, celltype) tuples
        "depth_count" : dict {celltype: {depth: count}}
        "leaf_order"  : list of leaf names
        "valid_types" : list of valid (non-NaN) cell type labels
    """
    from matplotlib import gridspec

    # Taking care of the palette
    palette = get_color_palette(clone_tdata, color_key=color, palette=palette)

    add_obs_to_tree(clone_tdata, keys=[color, "germ_layer"])
    py.pp.add_depth(clone_tdata)
    py.tl.ancestral_states(clone_tdata, keys=color, tree=clone_key, method="mode")
    flag_mixed_nodes(clone_tdata, clone_key, cell_type_key=color)
    add_node_attrs_to_edges(clone_tdata, keys=[color, "germ_layer"])

    # LCA depth computation — computed once, returned for reuse downstream
    points, depth_count, leaf_order, valid_types = compute_lca_depth_points(
        clone_tdata, clone_key, color, n_leaf=n_leaf
    )

    # Create a tree panel and a right-side point panel
    fig = plt.figure(figsize=(16, 40))
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 2], wspace=0.08)
    ax_tree = fig.add_subplot(gs[0])
    ax_points = fig.add_subplot(gs[1])

    pos = ax_points.get_position()
    new_bottom = 0.121
    new_height = 0.748
    ax_points.set_position([pos.x0, new_bottom, pos.width, new_height])

    # Plotting Tree
    py.pl.branches(
        clone_tdata,
        depth_key="depth",
        color=color,
        palette=palette,
        extend_branches=True,
        angled_branches=True,
        tree=clone_key,
        legend=False,
        ax=ax_tree,
    )
    py.pl.nodes(
        clone_tdata,
        color=color,
        palette=palette,
        nodes="leaves",
        size=10,
        legend=False,
        tree=clone_key,
        ax=ax_tree,
    )
    ax_tree.set_title(f"Lineage tree — {clone_key}", fontsize=14)

    # Dotplot
    plot_depth_dotplot(ax_points, points, valid_types, palette, leaf_order)

    # Histogram overlay
    plot_depth_histogram(ax_points, depth_count, valid_types, palette)

    plt.tight_layout()
    if output_path_plot is None:
        plt.show()
        result_ax = ax_tree
    else:
        plt.savefig(output_path_plot + clone_key + "_tree_with_scatter_histogram.svg", bbox_inches='tight')
        result_ax = None
    plt.close(fig)

    return {
        "ax_tree": result_ax,
        "points": points,
        "depth_count": depth_count,
        "leaf_order": leaf_order,
        "valid_types": valid_types,
    }

def compute_lca_depth_points(clone_tdata, clone_key, color, n_leaf=5):
    """
    Computes sliding-window LCA depth stats for each cell type.

    Returns
    -------
    points : list of (depth, leaf_position, celltype) tuples
    depth_count : dict {celltype: {depth: count}}
    leaf_order : list of leaf names
    valid_types : list of valid (non-NaN) cell type labels
    """
    py.tl.tree_distance(clone_tdata, tree=clone_key, metric="lca", key_added="lca")
    lca_depths = clone_tdata.obsp["lca_distances"]
    if hasattr(lca_depths, "toarray"):
        lca_depths = lca_depths.toarray()

    leaf_order = py.get.leaves(clone_tdata, tree=clone_key)
    n = len(leaf_order)

    # Precompute leaf -> row/col index once
    leaf_idx_map = {name: i for i, name in enumerate(clone_tdata.obs_names)}
    leaf_indices = np.array([leaf_idx_map[leaf] for leaf in leaf_order])

    # Reorder the full LCA matrix to leaf_order once
    lca_leaf = lca_depths[np.ix_(leaf_indices, leaf_indices)]

    cell_types = clone_tdata.obs.loc[leaf_order, color].astype(str).to_numpy()
    valid_types = [ct for ct in np.unique(cell_types) if ct not in {"nan", "None", "", "<NA>"}]

    # Assign each leaf a window id, build one small dataframe to group by (window, cell_type)
    window_id = np.arange(n) // n_leaf
    df = pd.DataFrame({"window_id": window_id, "cell_type": cell_types, "pos": np.arange(n)})
    df = df[df["cell_type"].isin(valid_types)]

    # Only windows with >=2 total leaves are valid (matches original's window-level skip)
    window_sizes = pd.Series(window_id).value_counts()
    valid_windows = set(window_sizes[window_sizes >= 2].index)

    points = []
    depth_count = {ct: {} for ct in valid_types}

    # groupby only visits (window, cell_type) pairs that actually occur —
    # no wasted iterations over empty combinations
    for (w, ct), group in df.groupby(["window_id", "cell_type"], sort=False):
        if w not in valid_windows or len(group) < 2:
            continue
        idx = group["pos"].to_numpy()
        middle_pos = min(w * n_leaf + 2, n - 1)

        sub = lca_leaf[np.ix_(idx, idx)]
        tri = np.triu_indices(len(idx), k=1)
        pair_depths = sub[tri]
        ancestor_depth = float(pair_depths.min()) if pair_depths.size else 0.0

        points.append((ancestor_depth, middle_pos, ct))
        depth_count[ct][ancestor_depth] = depth_count[ct].get(ancestor_depth, 0) + 1

    return points, depth_count, leaf_order, valid_types


def plot_depth_dotplot(ax_points, points, valid_types, palette, leaf_order):
    """
    Scatter plot of sliding-window ancestor depth vs leaf order, colored by cell type.
    """
    for ct in valid_types:
        subset = [p for p in points if p[2] == ct]
        if not subset:
            continue
        xs = [p[0] for p in subset]
        ys = [p[1] for p in subset]
        ax_points.scatter(xs, ys, color=palette.get(ct, "black"), s=40, alpha=0.8, label=ct)

    x_max = max([p[0] for p in points], default=1.0)
    ax_points.set_xlim(0, x_max + 1)
    ax_points.set_ylim(0, len(leaf_order) - 1)
    ax_points.set_xlabel("Ancestor depth")
    ax_points.set_ylabel("Leaf order (bottom -> top)")
    ax_points.set_title("Sliding-window ancestor-depth", fontsize=14)
    ax_points.grid(alpha=0.2)
    ax_points.legend(loc="lower left", frameon=False)

    return ax_points


def plot_depth_histogram(ax_points, depth_count, valid_types, palette):
    """
    Overlays a step histogram of ancestor-depth counts per cell type on a
    twin y-axis of ax_points. Returns the twin axis (ax_hist).
    """
    all_counts = [count for ct in valid_types for count in depth_count[ct].values()]
    max_count = max(all_counts, default=1)

    ax_hist = ax_points.twinx()

    for ct in valid_types:
        if not depth_count.get(ct):
            continue
        depths = list(depth_count[ct].keys())
        counts = list(depth_count[ct].values())
        if len(depths) == 1:
            bins = np.arange(min(depths) - 0.5, max(depths) + 1.5, 1)
        else:
            bins = np.arange(0, max(depths) + 2, 1)
        ax_hist.hist(
            depths,
            bins=bins,
            weights=counts,
            histtype='step',
            color=palette.get(ct, "black"),
            edgecolor=palette.get(ct, "black"),
            label=f"{ct} counts",
            stacked=True
        )

    ax_hist.set_ylim(0, max_count * 10)
    ax_hist.set_ylabel("Histogram share (scaled)")
    ax_hist.set_yticks([])

    return ax_hist

def get_color_palette(adata, color_key, palette=None, default_cmap="tab10"):
    """
    Build a {category: color} dict for adata.obs[color_key], staying consistent
    with any existing palette saved in adata.uns, and always including gray
    for NaN / missing values.

    Parameters
    ----------
    adata
        AnnData / TreeData object
    color_key : str
        obs column to build colors for
    palette : str | list | dict | None
        - None: look for adata.uns[f"{color_key}_colors"] first (scanpy convention);
          if not found, generate from `default_cmap`
        - str: either a key in adata.uns (a saved palette) or a matplotlib
          colormap name (e.g. "tab10", "tab20", "Set2")
        - list: ordered colors, zipped with sorted categories
        - dict: explicit mapping; any categories it's missing get auto-filled

    Returns
    -------
    dict mapping category -> color, plus "mixed" and NaN entries
    """
    categories = adata.obs[color_key].astype("category").cat.categories.tolist()
    resolved = {}

    if isinstance(palette, dict):
        resolved = dict(palette)
    elif isinstance(palette, list):
        resolved = dict(zip(categories, palette))
    elif isinstance(palette, str):
        if palette in adata.uns:
            resolved = dict(zip(categories, adata.uns[palette]))
        else:
            cmap = plt.get_cmap(palette)
            resolved = {cat: cmap(i % cmap.N) for i, cat in enumerate(categories)}
    elif palette is None:
        uns_key = f"{color_key}_colors"
        if uns_key in adata.uns:
            resolved = dict(zip(categories, adata.uns[uns_key]))
        else:
            cmap = plt.get_cmap(default_cmap)
            resolved = {cat: cmap(i % cmap.N) for i, cat in enumerate(categories)}

    # fill any categories still missing a color (e.g. dict/list didn't cover all)
    cmap = plt.get_cmap(default_cmap)
    for i, cat in enumerate(categories):
        if cat not in resolved:
            resolved[cat] = cmap(i % cmap.N)

    # reserved, always-consistent entries
    resolved["mixed"] = (0.6, 0.6, 0.6, 1.0)
    resolved[np.nan] = "lightgray"
    resolved["nan"] = "lightgray"   # covers the common case where obs got cast to str

    return resolved