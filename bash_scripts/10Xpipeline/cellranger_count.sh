#!/bin/bash
annotation=$1
module load bcl2fastq2/2.19.1
cellranger count  --id=$annotation \
                  --transcriptome=/data/neurogen/referenceGenome/Homo_sapiens/cellranger_references/hg38/transcriptome/refdata-gex-GRCh38-2020-A_Ensembl98_DL-2-25-2021/ \
                  --fastqs=../cellranger_fastqs/outs/fastq_path/ \
                  --sample=$annotation \
                  --include-introns \
                  --localcores=8 \
                  --localmem=64


