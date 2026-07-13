import scanpy as sc
import pandas as pd


def merge_robin_pijuan_atlas(
    robin_atlas, 
    pijuan_atlas,
    time, 
    celltype_label, 
    robin_time_label = "TimeStampMerged", 
    robin_celltype_label = "Annotation",
    pijuan_time_label = "stage",
    pijuan_celltype_label = "celltype_PijuanSala2019"
):
    '''
    Merging robin's celltype annotation into pijuan's old 2019 atlas, replacing only the cardiac celltype labels 
    and keeping every other type of label the same for background
    '''



    # subset pijuan to desired timepoint 
    pijuan_new_atlas = pijuan_atlas[pijuan_atlas.obs[pijuan_time_label] == time] 

    # Taking only the cells from 2019
    pijuan_atlas = pijuan_atlas[pijuan_atlas.obs[pijuan_celltype_label] != "New cells"]
    # setting index
    pijuan_atlas.obs.set_index('cell', inplace=True)

    # Subset Robin to pijuan cells only
    robin_atlas = robin_atlas[robin_atlas.obs["Dataset"] == "Pijuan.2019"]

    if robin_atlas.obs.index.isin(pijuan_atlas.obs.index).all() == False:
        raise Exception("Index values don't match. Make sure you set the cell index for both dataset")

    pijuan_atlas.obs[celltype_label] = robin_atlas.obs[robin_celltype_label].combine_first(pijuan_atlas.obs[pijuan_celltype_label])
    
    pijuan_atlas.layers['raw_counts'] = pijuan_atlas.X
    return pijuan_atlas


