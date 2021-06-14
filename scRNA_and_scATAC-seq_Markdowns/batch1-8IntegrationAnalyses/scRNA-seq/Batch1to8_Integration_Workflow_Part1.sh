#!/bin/bash
#SBATCH -t 96:00:00
#SBATCH --mem=200G
#SBATCH -p medium
#SBATCH -e logs/Rscript-%J.log
conda activate r4-base
Rscript Batch1to8_Integration_Workflow_Part1.R
