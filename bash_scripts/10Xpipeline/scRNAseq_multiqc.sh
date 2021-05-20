#!/bin/bash
module load fastqc/0.11.8
module load multiqc/1.7
#reading in sample name passed from submission script
annotation=$1
#setting up tempfiles and outfiles for files consisting of all merged fastqs for given sample
mergedI1file=${annotation}.index1.fastq
zipmergedI1file=${annotation}.index1.fastq.gz
mergedI2file=${annotation}.index2.fastq
zipmergedI2file=${annotation}.index2.fastq.gz
#setting up tempfiles and outfiles for files consisting exclusively of merged read fastqs (R1 & R2, index fastqs - I1 & I2 - are excluded) for given sample
mergedR1file=${annotation}.read1.fastq
zipmergedR1file=${annotation}.read1.fastq.gz
mergedR2file=${annotation}.read2.fastq
zipmergedR2file=${annotation}.read2.fastq.gz
#specifying the report filename
reportname=${annotation}_multiqc_report.html
#naming the outdir
outdir=${annotation}_fastqc
#merging all index 1 (forward index) files for sample and zipping output
zcat ../cellranger_fastqs/outs/fastq_path/*/${annotation}_*_*_I1_* > $mergedI1file
gzip $mergedI1file
#merging all index 2 (reverse index) files for sample and zipping output
zcat ../cellranger_fastqs/outs/fastq_path/*/${annotation}_*_*_I2_* > $mergedI2file
gzip $mergedI2file
#merging all forward read (excluding index1) files for sample and zipping output
zcat ../cellranger_fastqs/outs/fastq_path/*/${annotation}_*_*_R1_* > $mergedR1file
gzip $mergedR1file
#merging all reverse read (excluding index2) files for sample and zipping output
zcat ../cellranger_fastqs/outs/fastq_path/*/${annotation}_*_*_R2_* > $mergedR2file
gzip $mergedR2file
#making the outdir
mkdir $outdir
#running fastqc on the merged forward and reverse read and index files
fastqc --outdir $outdir $zipmergedR1file $zipmergedR2file $zipmergedI1file $zipmergedI2file ../cellranger_counts/${annotation}/outs/possorted_genome_bam.bam
multiqc --outdir multiQC_reports.dir --filename $reportname $outdir
#removing merged files
rm $zipmergedI1file
rm $zipmergedI2file
rm $zipmergedR1file
rm $zipmergedR2file 


