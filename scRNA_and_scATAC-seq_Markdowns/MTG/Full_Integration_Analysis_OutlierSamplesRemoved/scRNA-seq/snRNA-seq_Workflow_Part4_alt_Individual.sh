#!/bin/bash
#SBATCH -t 96:00:00
#SBATCH --mem=900G
#SBATCH -p highmem
module load conda2/4.2.13
source activate /home/jap0606/.conda/envs/condamamba/envs/r4-base
rds_path=$(head -n 1 Files/SubclusterIdents.txt)
cluster=GLU_Neurons_1
Rscript snRNA-seq_Workflow_Part4_alt_Individual.R --SeuratObject $rds_path --ClusterIdent $cluster
