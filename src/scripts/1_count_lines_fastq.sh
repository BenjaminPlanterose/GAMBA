#!/bin/bash

FILES=$(ls *.fastq.gz)

for i in $FILES
do
	echo $i >> counts.txt
	echo $(zcat $i | wc -l)/4|bc >> counts.txt
done
