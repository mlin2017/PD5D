#!/bin/bash
#SBATCH -t 96:00:00
#SBATCH --mem=150G
#SBATCH -p medium
#SBATCH -e logs/Rscript-%J.log
Rscript Batch1to8_Integration_Workflow_Part3.R
