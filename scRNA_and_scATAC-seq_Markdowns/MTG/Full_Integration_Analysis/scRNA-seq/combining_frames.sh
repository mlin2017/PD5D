#!/bin/bash
#SBATCH -t 96:00:00
#SBATCH --mem=700G
#SBATCH -p highmem
#SBATCH -e logs/Combinescript-%J.log
module load conda2/4.2.13
source activate /home/jap0606/.conda/envs/condamamba/envs/cgatflow
for i in $(ls /n/scratch3/users/j/jap0606/FullIntegration/batch*/matrix.mtx); do
    tail -n +4 $i >> /n/scratch3/users/j/jap0606/FullIntegration/PreCombinedFiles/matrix.mtx
done 
gzip /n/scratch3/users/j/jap0606/FullIntegration/PreCombinedFiles/matrix.mtx

#paste --delimiters=' ' /n/scratch3/users/j/jap0606/FullIntegration/batch*/matrix.mtx | tail -n +3 >> /n/scratch3/users/j/jap0606/FullIntegration/PreCombinedFiles/matrix.mtx
#gzip /n/scratch3/users/j/jap0606/FullIntegration/PreCombinedFiles/matrix.mtx
