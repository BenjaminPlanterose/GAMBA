#!/bash/bin

# Manual page
if [ "$1" == "-h" ] || [ "$1" == "--h" ] || [ "$1" == "--help" ] || [ "$1" == "-help" ]; then
  echo "Usage: $0 -f DIR -p FILE -g FILE -t DIR [-q INT -S INT -L INT -c INT]" >&2
  echo
  echo "   -f           Directory containing paired-end fastq files."
  echo "   -p           Primer sequences in FASTA format."
  echo "   -g           Directory containing genome/Amplicon sequences in FASTA/Multi-FASTA format and its bisulfite converted counterpart (see bismark_genome_preparation)."
  echo "   -t           Path the trimmomatic jar executable."
  echo "   -q           Read quality filter. Default is 20 (i.e. Q20)."
  echo "   -S           Low-end size filter. Default is 0 bp."
  echo "   -L           High-end size filter. Default is 1000 bp."
  echo "   -c           Number of cores to employ. Default is 1 core."
  echo  
  exit 0
fi

# Parse arguments
while getopts ":f:p:g:t:q:S:L:c:" opt; do
  case $opt in
    f) fastq_path="$OPTARG"
    ;;
    p) primers="$OPTARG"
    ;;
    g) genome="$OPTARG"
    ;;
    t) trimmomatic_path="$OPTARG"
    ;;
    q) QC="$OPTARG"
    ;;
    S) size_S="$OPTARG"
    ;;
    L) size_L="$OPTARG"
    ;;
    c) nCores="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac

  case $OPTARG in
    -*) echo "Option $opt needs a valid argument"
    exit 1
    ;;
  esac
done

# Set parameters with default values if user has not input them
if [ -z "$nCores" ]; then
  nCores=1
fi

if [ -z "$size_L" ]; then
  size_L=1000
fi

if [ -z "$size_S" ]; then
  size_S=0
fi

if [ -z "$QC" ]; then
  QC=20
fi


# Make sure all mandatory arguments have been input
if [ -z "$fastq_path" ]; then
  printf "***********************************\n"
      printf "* Error: Directory to fastq files not supplied (tag -f).*\n"
      printf "***********************************\n"
      exit 1
fi

if [ -z "$primers" ]; then
  printf "**********************************\n"
      printf "* Error: Primers FASTA file not supplied (tag -p).*\n"
      printf "**********************************\n"
      exit 1
fi

if [ -z "$genome" ]; then
  printf "**********************************\n"
      printf "* Error: Genome FASTA file not supplied (tag -g).*\n"
      printf "**********************************\n"
      exit 1
fi

if [ -z "$trimmomatic_path" ]; then
  printf "**********************************\n"
      printf "* Error: Trimmomatic path not supplied (tag -t).*\n"
      printf "**********************************\n"
      exit 1
fi

# Find path to GAMBA
MY_PATH=$(dirname "$0")
#MY_PATH=$(pwd)

# Print selected arguments
printf "Argument fastq_path is %s\n" "$fastq_path"
printf "Argument primers is %s\n" "$primers"
printf "Argument genome is %s\n" "$genome"
printf "Argument trimmomatic_path is %s\n" "$trimmomatic_path"
printf "Argument QC is %s\n" "$QC"
printf "Argument size_S is %s\n" "$size_S"
printf "Argument size_L is %s\n" "$size_L"
printf "Argument nCores is %s\n" "$nCores"
printf "Argument GAMBA_PATH is %s\n" "$MY_PATH"

# Begin GAMBA pipeline
cd $fastq_path
bash $MY_PATH/scripts/0_main.sh $fastq_path $primers $genome $trimmomatic_path $QC $size_S $size_L $nCores $MY_PATH


