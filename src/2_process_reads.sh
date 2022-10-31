#!/bash/bin


TRIM=$1
primer=$2
fastq=$3
QC=$4
m=$5
M=$6
nThread=$7

R1=($( ls $fastq/*R1*.fastq.gz -1v))
R2=($( ls $fastq/*R2*.fastq.gz -1v))
l=${#R1[@]}

echo ${R1[@]}
echo ${R2[@]}

mkdir trimmomatic
cd trimmomatic
mkdir cutadapt

for ((i=0;i<=$l-1;i+=1))
do
 preindex=$(echo ${R1[$i]} | rev | cut -d'/' -f 1 | rev)
 index="$(cut -d'_' -f2 <<< ${preindex})" 
 echo ${R1[$i]} ${R2[$i]} $index
 java -jar $TRIM/trimmomatic-0.39.jar PE ${R1[$i]} ${R2[$i]} -summary summary_${index}.txt -threads $nThread -baseout cut_${index}.fastq.gz ILLUMINACLIP:$primer/primers.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
 cutadapt -m $m -M $M -q $QC -j $nThread -o cutadapt/trim_cut_${index}_R1.fastq.gz -p cutadapt/trim_cut_${index}_R2.fastq.gz cut_${index}_1P.fastq.gz cut_${index}_2P.fastq.gz > cutadapt/cutadapt.log${index}
done



