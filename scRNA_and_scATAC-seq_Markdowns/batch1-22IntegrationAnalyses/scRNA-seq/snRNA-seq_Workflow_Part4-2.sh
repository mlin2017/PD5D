#!/bin/bash
#SBATCH -t 96:00:00
#SBATCH --mem=200G
#SBATCH -p medium
module load conda2/4.2.13
source activate /home/jap0606/.conda/envs/condamamba/envs/r4-base
rds_path=$(head -n 1 Files/ClusterIdents.txt)
cluster=$clusterident
Rscript snRNA-seq_Workflow_Part4-2.R --SeuratObject $rds_path --ClusterIdent $cluster
