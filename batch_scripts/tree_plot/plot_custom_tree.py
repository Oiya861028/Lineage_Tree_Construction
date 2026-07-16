import scanpy as sc
import treedata as td
import pycea as py
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np
import pandas as pd
import gc
import sys
sys.path.append("/project/imoskowitz/yubin/Lineage_Tree_Construction")
from lineage_utilities.tree_functions import *

data_dir = "output_data"
plot_dir = "output_plot"
base_path = "/project/imoskowitz/yubin/Lineage_Tree_Construction/"
output_path_data = base_path+data_dir+"/"
output_path_plot = base_path+plot_dir+"/Trees/Robin_Pijuan/Manual_Cardiac_Annotation/E8_5/"
adata_fname = "Processed_data/E8_5.h5td"


adata = td.read_h5td(output_path_data+adata_fname)
adata.X = adata.layers["raw_counts"]

adata.obs['cardiac_labels'] = adata.obs['cardiac_labels'].astype('category')
adata.obs['tree'] = adata.obs['tree'].astype('category')
adata.obs['germ_layer'] = adata.obs['germ_layer'].astype('category')

tree_sizes = adata.obs['tree'].value_counts().sort_values(ascending=False)
print(tree_sizes.head(10))

# flag trees likely to be memory-heavy (dense O(n^2) LCA matrix)
LARGE_TREE_THRESHOLD = 2000  # adjust based on what crashed before
large_trees = tree_sizes[tree_sizes > LARGE_TREE_THRESHOLD].index.tolist()
print("Large trees (may need special handling):", large_trees)



failed_trees = []

output_csv = output_path_plot + "depth_points_all_trees.csv"
csv_path = Path(output_csv)
if csv_path.exists():
    csv_path.unlink()  # start fresh; remove if you want to append across runs

for clone_key in adata.obs['tree'].unique():
    if pd.isna(clone_key) or clone_key == 'nan': # To catch cells that do not have a tree
        continue

    n_cells = tree_sizes.get(clone_key, 0)
    print(f"starting {clone_key} ({n_cells} cells)")
    
    try: 
        clone_tdata = adata[adata.obs["tree"] == clone_key].copy()

        # plotting + LCA computation happen once, together
        result = plot_tree_custom(clone_tdata, clone_key, color="cardiac_labels", output_path_plot=output_path_plot)

        points = result["points"]
        leaf_order = result["leaf_order"]

        rows = [{"tree": clone_key, "cell_type": ct, "depth": d, "pos": p, "n_leaves": len(leaf_order)}
                for d, p, ct in points]
        df_chunk = pd.DataFrame(rows)
        df_chunk.to_csv(output_csv, mode='a', header=not csv_path.exists(), index=False)
        del clone_tdata, result, points
    except MemoryError as e:
        print(f"MEMORY ERROR on {clone_key}, skipping: {e}")
        failed_trees.append(clone_key)
    except Exception as e:
        print(f"ERROR on {clone_key}, skipping: {e}")
        failed_trees.append(clone_key)
    finally:
        gc.collect()

    print(f"complete {clone_key}")

plt.close('all')  # belt-and-suspenders in case anything else left a figure open
print("Failed trees:", failed_trees)