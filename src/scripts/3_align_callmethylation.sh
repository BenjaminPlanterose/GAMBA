#!bin/bash


genome=$1
fastq=$2
nThread=$3


R1=($( ls $fastq/*R1*.fastq.gz -1v))
R2=($( ls $fastq/*R2*.fastq.gz -1v))
l=${#R1[@]}


for ((i=0;i<=$l-1;i+=1))
do
 preindex=$(echo ${R1[$i]} | rev | cut -d'/' -f 1 | rev)
 index="$(cut -d'_' -f3 <<< ${preindex})" 
 echo ${R1[$i]} ${R2[$i]} $index
 bismark --genome $genome --non_directional --multicore $nThread -D 15 -R 2 -1 $fastq/trim_cut_${index}_R1.fastq.gz -2 $fastq/trim_cut_${index}_R2.fastq.gz > log${index}.txt
done

bismark2report
bismark2summary

bams=($( ls *.bam -1v))
l=${#bams[@]}

cd ..
cd 3_methylation

for ((i=0;i<=$l-1;i+=1))
do
 preindex=$(echo ${bams[$i]} | rev | cut -d'/' -f 1 | rev)
 index="$(cut -d'_' -f3 <<< ${preindex})" 
 echo ${bams[$i]} $index
 bismark_methylation_extractor -p --cutoff 1 --gzip --multicore $nThread --comprehensive --bedGraph --CX_context --zero_based --cytosine_report --no_overlap --genome_folder $genome ../2_alignment/${bams[$i]} > log${index}.txt
done

mkdir bedGraphs
mv *cov.gz bedGraphs






