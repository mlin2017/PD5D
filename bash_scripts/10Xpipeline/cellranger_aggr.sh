#/bin/bash
module load python/3.8.2
python mkaggrtable.py 
cellranger aggr  --id=Batch8_050521 \
                 --csv=aggrtable.csv \
                 --normalize=none


