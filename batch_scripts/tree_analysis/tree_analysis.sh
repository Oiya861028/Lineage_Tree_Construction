#!/bin/bash
#SBATCH --job-name=tree_analysis_exclude_top3_largest_tree
#SBATCH --account=pi-imoskowitz
#SBATCH --output=/project/imoskowitz/yubin/Lineage_Tree_Construction/batch_scripts/tree_analysis/E8_5_tree_run_2_only_excluded_from_run_1%j.out 
#SBATCH --error=/project/imoskowitz/yubin/Lineage_Tree_Construction/batch_scripts/tree_analysis/E8_5_tree_%j.err 
#SBATCH --partition=bigmem
#SBATCH --mem=256G
#SBATCH --cpus-per-task=16
#SBATCH --time=36:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yubin@rcc.uchicago.edu


module load python
source activate Lineage_Tracing
# module unload python

python /project/imoskowitz/yubin/Lineage_Tree_Construction/batch_scripts/tree_analysis/run_tree_analysis.py
