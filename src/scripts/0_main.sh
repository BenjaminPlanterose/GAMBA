#!/bin/bash

#################################################### Params

where=$1
PRIMER=$2
genome=$3
TRIM=$4
QC=$5
m=$6
M=$7
nThread=$8
scripts="$9/scripts"

#################################################### Main

cd $where
mkdir 0_FASTQ
mv *.fastq.gz ./0_FASTQ

# 1. Run QC
cd ./0_FASTQ/
bash $scripts/1_fastqc.sh $nThread
bash $scripts/1_count_lines_fastq.sh

# 2. Filtering
cd ..
mkdir 1_trimmed_reads
cd 1_trimmed_reads
bash $scripts/2_process_reads.sh $TRIM $PRIMER $where/0_FASTQ $QC $m $M $nThread
where2="$where/1_trimmed_reads/trimmomatic/cutadapt"

# 3. Run QC post-cuttin/trimming
cd $where2
bash $scripts/1_fastqc.sh $nThread
bash $scripts/1_count_lines_fastq.sh

# 4. Align and perform DNA methylation calling
cd $where
mkdir 2_alignment
mkdir 3_methylation
cd 2_alignment
bash $scripts/3_align_callmethylation.sh $genome $where2 $nThread
Rscript $scripts/4_methylation_analysis.R $where/3_methylation/bedGraphs/ $nThread

# 5. Sort and index bam files for visualization in IGV
where3="$where/2_alignment"
bash $scripts/5_index_bam_visual.sh $where3
bash $scripts/5_count_alignment_bam.sh
