#!/bin/bash
#SBATCH -t 96:00:00
#SBATCH --mem=250G
#SBATCH -p medium
#SBATCH -e logs/Rscript-%J.log
module load conda2/4.2.13
source activate /home/jap0606/.conda/envs/condamamba/envs/r4-base
Rscript making_combined_matrix.R
