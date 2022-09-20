#!/bin/bash
#SBATCH -t 96:00:00
#SBATCH --mem=50G
#SBATCH -p medium
#SBATCH -e logs/Rscript-%J.log
for Ident in $(tail -n +2 Files/ClusterIdents.txt); do
    sbatch --export=clusterident=$Ident DE_Test.sh
done

