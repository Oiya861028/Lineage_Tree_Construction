import scanpy as sc
import pandas as pd


def merge_robin_pijuan_atlas(
    robin_atlas, 
    pijuan_atlas,
    time, 
    celltype_label, 
    robin_celltype_label = "Annotation",
    pijuan_time_label = "stage",
    pijuan_celltype_label = "celltype_PijuanSala2019"
):
    '''
    Merging robin's celltype annotation into pijuan's old 2019 atlas, replacing only the cardiac celltype labels 
    and keeping every other type of label the same for background

    Parameters 
    ----------
    robin_atlas
        Robin's time specific adata containing celltype annotation

    pijuan_atlas
        Pijuan's full atlas with all time, will be subsetted with `time` parameter
    
    celltype_label
        the obs label that will store the merged celltype, should be the same name that you would pass into celltypist training
    
    robin_celltype_label
        Obs name for robin's celltype atlas
    
    pijuan_time_label
        Obs name for Pijuan's cell time
    
    pijuan_celltype_label 
        Obs name for Pijuan's old dataset celltype label

    Return
    ------
    Adata containing robin's cardiac labeling with only pijuan's old cells 
    '''



    # subset pijuan to desired timepoint 
    pijuan_atlas = pijuan_atlas[pijuan_atlas.obs[pijuan_time_label] == time] 

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


