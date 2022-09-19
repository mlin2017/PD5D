#!/bin/bash
#SBATCH -t 96:00:00
#SBATCH --mem=200G
#SBATCH -p medium
module load conda2/4.2.13
source activate /home/jap0606/.conda/envs/condamamba/envs/r4-base
rds_path=$(head -n 1 Files/ClusterIdents.txt)
cluster=$clusterident
Rscript Batch1to8_Integration_Workflow_Part4.R --SeuratObject $rds_path --ClusterIdent $cluster
