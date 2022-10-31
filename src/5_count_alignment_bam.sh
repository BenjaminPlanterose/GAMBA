#!/bin/bash


FILES=$(ls *sorted.bam)

for i in $FILES
do
	echo $i >> counts.txt
	samtools view $i -c >> counts.txt
done

