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
bash /home/ultron/Git/GAMBA/src/GAMBA.sh -f /home/ultron/Git/GAMBA/example/fastqs \
-p /home/ultron/Git/GAMBA/example/primers/primers.fa -g /home/ultron/Git/GAMBA/example/genome/ \
-t /home/ultron/opt/Trimmomatic-0.39/trimmomatic-0.39.jar -q 20 -S 100 -L 230 -c 1
```