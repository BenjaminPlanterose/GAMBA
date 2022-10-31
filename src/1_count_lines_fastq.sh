#!/bin/bash

FILES=$(ls *.fastq.gz)

for i in $FILES
do
	echo $i >> counts.txt
	#zcat $i | grep -c '^@' >> counts.txt
	echo $(zcat $i | wc -l)/4|bc >> counts.txt
done

