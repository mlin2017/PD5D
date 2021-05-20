#!/bin/bash
for i in $(cat annotations*.csv | awk -F "," '(NR>1){print $2}'); do
    bsub -q big -e ${i}_count.log "sh cellranger_count.sh $i"
done
