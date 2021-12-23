#!/bin/bash
for filename in $(ls *.fastq.gz);do
    sample=$(echo $filename | grep -o -E BN[0-9]+\|H[0-9]+\|BRI[0-9]+)
    disease_case=$(grep $sample PD5D_Sequencing_Tracking_Sheets_scRNA_seq_QC_filled.tsv | awk '{print $4}')
    tissue=$(grep $sample PD5D_Sequencing_Tracking_Sheets_scRNA_seq_QC_filled.tsv | awk '{print $3}')
    batch="batch"$(grep $sample PD5D_Sequencing_Tracking_Sheets_scRNA_seq_QC_filled.tsv | awk '{print $11}')
    suffix=$(echo $filename | grep -o -E L.*_R[1-2])
    new_filename=${disease_case}"_"${sample}"_"${tissue}"_"${batch}_${suffix}.fastq.gz
    rename $filename $new_filename * 
done
