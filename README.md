![Alt text](img/GAMBA.png?raw=true "logo")

# Genomic Analysis by MPS on Bisulfite-converted Amplicons (GAMBA)

## Requirements

```
Operating system: tested on Ubuntu 18.04.6 LTS (Bionic Beaver)
R: tested on R version 4.2.1 (2022-06-23) -- "Funny-Looking Kid"
```

## Dependencies

* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [MultiQC](https://multiqc.info/)
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
* [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)
* [samtools](http://www.htslib.org/)
* R-packages
	* [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
	* [R.utils](https://cran.r-project.org/web/packages/R.utils/index.html)

## Example


```bash
cd <path_to_GAMBA>/example/fastqs/
bash <path_to_GAMBA>/src/GAMBA.sh -f <path_to_GAMBA>/example/fastqs \
-p <path_to_GAMBA>/example/primers/primers.fa -g <path_to_GAMBA>/example/genome/ \
-t <path_to_trimmomatic_jar_executable> -q 20 -S 100 -L 230 -c 1
```