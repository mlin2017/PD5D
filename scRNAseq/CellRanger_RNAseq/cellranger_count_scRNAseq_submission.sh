#!/bin/bash
for i in $(cat annotations_021721.csv | awk -F "," '(NR>1){print $2}'); do
    bsub -q big -e {$i}_cellrangercounts.log "sh cellranger_counts.sh $i"
done
