#!/bin/bash
#SBATCH -t 96:00:00
#SBATCH --mem=700G
#SBATCH -p highmem
#SBATCH -e logs/pythoncombinescript-%J.log
module load conda2/4.2.13
source activate /home/jap0606/.conda/envs/condamamba/envs/cgatflow
python combining_frames.py 
