#!/bin/bash
#SBATCH --job-name=tree_plot
#SBATCH --account=pi-imoskowitz
#SBATCH --output=/project/imoskowitz/yubin/Lineage_Tree_Construction/batch_scripts/tree_plot/E8_5_treeplot_%j.out 
#SBATCH --error=/project/imoskowitz/yubin/Lineage_Tree_Construction/batch_scripts/tree_plot/E8_5_treeplot_%j.err 
#SBATCH --partition=bigmem
#SBATCH --mem=256G
#SBATCH --cpus-per-task=16
#SBATCH --time=36:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yubin@rcc.uchicago.edu


module load python
source activate Lineage_Tracing
# module unload python

python /project/imoskowitz/yubin/Lineage_Tree_Construction/batch_scripts/tree_plot/plot_custom_tree.py
