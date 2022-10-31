#!/bash/bin

nThread=$1

# Run fastqc
mkdir QC
fastqc ./*.fastq.gz -o ./QC/ -t $nThread
cd QC
multiqc ./*fastqc.zip
