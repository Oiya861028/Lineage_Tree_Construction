#!/bin/bash
#SBATCH --job-name=robin_pijuan_E8.5
#SBATCH --account=pi-imoskowitz
#SBATCH --output=/project/imoskowitz/yubin/Lineage_Tree_Construction/output_data/Celltypist/celltypist_%j.out 
#SBATCH --error=/project/imoskowitz/yubin/Lineage_Tree_Construction/output_data/Celltypist/celltypist_%j.err
#SBATCH --partition=bigmem
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=36:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yubin@rcc.uchicago.edu


module load python
source activate Lineage_Tracing
module unload python

python /project/imoskowitz/yubin/Lineage_Tree_Construction/Moslin/celltypist_script/run_robin_pijuan_celltypist_training.py
