#!/bin/bash/

where=$1

bams=($( ls $where/*.bam -1v))
l=${#bams[@]}
#index=( $(seq 1 $l ) )

cd $where

for ((i=0;i<=$l-1;i+=1))
do
 preindex=$(echo ${bams[$i]} | rev | cut -d'/' -f 1 | rev)
 index="$(cut -d'_' -f3 <<< ${preindex})" 
 echo ${index} 
 samtools sort ${bams[$i]} > ${index}_sorted.bam
 samtools index ${index}_sorted.bam
done



