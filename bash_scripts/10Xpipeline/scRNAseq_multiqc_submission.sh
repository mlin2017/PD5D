#!/bin/bash
for i in $(cat annotations_*.csv | awk -F "," '(NR>1){print $2}'); do
    bsub -q big -e ${i}_error.log "sh scRNAseq_multiqc.sh $i"
done
