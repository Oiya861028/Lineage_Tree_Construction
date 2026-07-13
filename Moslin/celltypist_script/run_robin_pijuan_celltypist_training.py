import sys
import scanpy as sc

# Importing celltypist training function from Katie
sys.path.append("/project/imoskowitz/shared/sequencing.processed/Smo_null_snRNAseq2025/3-Brain_system/Brain_system_analysis_python") 
from src.III_celltype_annotation.annotate_celltypist import (
    train_celltypist_model
)
# Importing merging atlas function from Yubin 
sys.path.append("/project/imoskowitz/yubin/Lineage_Tree_Construction/Moslin/celltypist_script")
from merging_atlas import (
    merge_robin_pijuan_atlas # Merging Script designed on Jul 10 2026, some assumption was made for obs names from both atlas. 
)

output_path_model = "/project/imoskowitz/yubin/Lineage_Tree_Construction/output_data/Celltypist/Models/Robin_Pijuan/E8_5/" # Change this 
celltype_label = "merged_robin_pijuan_celltype"

# Import Robin atlas
robin_atlas = sc.read_h5ad(
    "/project/imoskowitz/kdreyer/lab_datasets/002_Cardio_mesodermal_atlas/Cardio-mesodermal_atlas_formatted_E85.h5ad" # and this 
)

# Import Pijuan full atlas (Will be subsetted in the merging function).
pijuan_atlas = sc.read_h5ad(
    "/project/imoskowitz/kdreyer/celltypist_models/source_data/gastrulation_extended_E75_E775_E80_E825_E85_E875.h5ad"
)


# robin_celltype_label name: E8_0 cell_type, E8_5 Annotation
merged_atlas = merge_robin_pijuan_atlas(robin_atlas, pijuan_atlas, time = "E8.5", robin_celltype_label= "Annotation", celltype_label = celltype_label) # and time 
merged_atlas_fname = "Robin_Pijuan_annotation_E8_5"

model, model_fname = train_celltypist_model(
    adata_atlas=merged_atlas, adata_atlas_fname=merged_atlas_fname,
    celltype_label=celltype_label, output_path_model=output_path_model,
    top_genes=1000
)

print(model_fname)