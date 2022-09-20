#!/bin/bash
#SBATCH -t 96:00:00
#SBATCH --mem=50G
#SBATCH -p medium
#SBATCH -e logs/Rscript-%J.log
for Ident in $(tail -n +2 Files/ClusterIdentsSubset.txt); do
    sbatch --export=clusterident=$Ident snRNA-seq_Workflow_Part4_alt_Individual.sh
done

