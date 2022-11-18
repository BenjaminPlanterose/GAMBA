![Alt text](img/GAMBA.png?raw=true "logo")

# Genomic Analysis by MPS on Bisulfite-converted Amplicons (GAMBA)

## Requirements

```
Operating system: tested on Ubuntu 18.04.6 LTS (Bionic Beaver)
R: tested on R version 4.2.1 (2022-06-23) -- "Funny-Looking Kid"
```

For issues arising while running GAMBA, either [report an issue](https://github.com/BenjaminPlanterose/GAMBA/issues) or simply contact b.planterosejimenez at erasmusmc.nl.

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

The PATH for dependencies *FastQC*, *MultiQC*, *cutadapt*, *bismark* and *samtools* are assumed to be exported in the [.bashrc](https://linuxhint.com/export-a-path-in-bashrc/) file.
The PATH for dependency *Trimmomatic* instead must be supplied as an input.

## Example

An example dataset is available at "<path_to_GAMBA>/example/". To run GAMBA on it, run the following bash commands:

```bash
cd <path_to_GAMBA>/example/fastqs/
bash <path_to_GAMBA>/src/GAMBA.sh -f <path_to_GAMBA>/example/fastqs \
-p <path_to_GAMBA>/example/primers/primers.fa -g <path_to_GAMBA>/example/genome/ \
-t <path_to_trimmomatic_jar_executable> -q 20 -S 100 -L 230 -c 1
```

The expected output is available as "<path_to_GAMBA>/example/expected_output.zip".


## Details for its implementation

### Sample names

GAMBA assumes certain formating on the name of the sequencing files. Specifically, the following name processing is performed: 
```bash
file_name="<PATH>/84_S84_L001_R1_001.fastq.gz"
preindex=$(echo $file_name | rev | cut -d'/' -f 1 | rev) # 84_S84_L001_R1_001.fastq.gz
index="$(cut -d'_' -f2 <<< ${preindex})" # S84
```

Thus, make sure that this convention is used on your specific samples.

### Genome preparation

The genome file (-g) is a multi-FASTA file where each sequence corresponds to a different amplicon. 
All nucleotide sequence must contain primer-binding sites; additionally, it is critical to add NN at the beginning and at the end of each sequence (i.e. NN pre- and post- padding).
If not, the Bismark bisulfite aligner (in principle, not designed for amplicon bisulfite sequencing but rather whole genome bisulfite sequencing) will drop out most of the alignment.
Additionally, it is necessary to perform genome preparation as per [Bismark user-guide](https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html):

```bash
samtools index <genome_file>
bismark_genome_preparation <genome_file>
```

### Primer preparation

The primer file (-p) is as well a multi-FASTA file. In our case, to prepare out file, we started from a "primers.txt" file with the following format:

| cg_name       | Primer                  | Direction |
| ------------- | ----------------------- | --------- |
|cg13039251	|TGTGTAAGTTAGTTTGTGTT	  | Fwd       |
|cg15693572	|TTAAAATGTAGATTAGGGGA	  | Fwd       |
|cg23576855	|TAGGGTTGTTTTTTTAGAG	  | Fwd       |
|cg01940273	|GATAAAGTTTGGTTTTTTGG	  | Fwd       |
|cg03636183	|TTTATTAGTAGTATGGTGGA	  | Fwd       |
|cg05575921	|ATTGTTTATTTTTGAGAGGG	  | Fwd       |
|cg12876356	|TGGGTATATTGATTTTTTT	  | Fwd       |
|cg21566642	|TAGTTGGGGTTTTTGTATTTAG	  | Fwd       |
|cg22132788	|TATATTGTTAGGGGTGAGT	  | Fwd       |
|cg05951221	|TTGGTTGTTAGGAGGTT	  | Fwd1      |
|cg05951221	|TTGGTTGTTAGGAGGTC	  | Fwd2      |
|cg06126421	|TTATGGTAATTGTTTTGGAG	  | Fwd       |
|cg12803068	|TTTTGTTGATAGGGGGAA	  | Fwd       |
|cg09935388	|TTAGTGAGAGGTTGTATTT	  | Fwd1      |
|cg09935388	|TTAGTGAGAGGTTGTATTC	  | Fwd2      |
|cg13039251	|AACTATCTCCCTATTTTCTA	  | Rv        |
|cg15693572	|CTCAACCACATTATCATAAAACA  | Rv        |
|cg23576855	|CCCTTCTTAATTACAATAAAC	  | Rv        |
|cg01940273	|ATTACATCTCTCTTCCCTT	  | Rv        |
|cg03636183	|ACCAAATCTATACCAATAAC	  | Rv        |
|cg05575921	|AACTCTATACCTCCAAAA	  | Rv        |
|cg12876356	|AATCTATTTACTATTCTACCTCC  | Rv        |
|cg21566642	|CTTAAATACTTAACCTCCT	  | Rv        |
|cg22132788	|CCAAAACAAAATAAAAAAAC	  | Rv        |
|cg05951221	|ACTTCTCTCAAAAAACAA	  | Rv        |
|cg06126421	|AATATTTCCCCTTTTATCCA	  | Rv        |
|cg12803068	|ACAATAACACATACAAAAACT	  | Rv        |
|cg09935388	|CCAACAAATATATCTAAAAACC	  | Rv        |

We generated a multi-FASTA file with the following R-script:

```r
library(data.table)
library(seqinr)
primers = as.data.frame(fread("primers.txt"))
primers_RC = primers
primers_RC$Primer = sapply(1:nrow(primers_RC), function(x) RC(primers_RC$Primer[x]))
primers$orient = "direct"
primers_RC$orient = "RC"
PRIMERS = rbind(primers, primers_RC)
LIST = lapply(1:nrow(PRIMERS), function(x) PRIMERS$Primer[x])
Names = paste(PRIMERS$`cg Name`, PRIMERS$Direction, primers$orient, sep = "_")
write.fasta(LIST, names = Names, file.out = "primers.fa")
```
