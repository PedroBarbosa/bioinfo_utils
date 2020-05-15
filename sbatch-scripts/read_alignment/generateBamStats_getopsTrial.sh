#!/bin/bash
usage="$(basename "$0") input_file mapper outdir [-h] [-s] [-o]. Script to generate an automated sbatch file to run a mapping stats analysis in a slurm-based cluster environment.
 where:
    input_file	Text file where each line represents the path for the SAM/BAM files to be processed.
    mapper	Aligner used to generate the alignments. Available options: Tophat2, Hisat2, STAR, bowtie2, bwa-mem.
    outdir      Directory to write the output file. 
    Optional arguments:
    -h  show this help text
    -o  name of the output file. Default:'mappingStats.txt'
    -s  flag to indicate that the reads are single-end. (default paired-end)"


while getopts ':hso:' option; do
  case "$option" in
    "h") echo "$usage"
       exit
       ;;
    "s") printf "Single end flag (-${option}) handling is not ready to be used for now. Exiting\n"
       echo "$usage"
       exit 1
       ;;
    "o") echo "-$option detected. Output file will be named $OPTARG"
         output_file=$OPTARG
       ;;
    ":")
       echo "Missing argument for -%s\n" "$OPTARG"
       echo "$usage"
       exit 1
       ;;
    "?") printf "illegal option: -%s\n" "$OPTARG" 
       echo "$usage" 
       exit 1
       ;;
  esac
#  echo "$OPTIND"
#  shift $((OPTIND - 1)) 

done

ARG1=${@:$OPTIND:1}
ARG2=${@:$OPTIND+1:1}
ARG3=${@:$OPTIND+1+1:1}
echo "first: $ARG1"
echo "second: $ARG2"
echo "third: $ARG3"
if [[ $# -ne 3 ]]; then
    printf "$0: ERROR Three arguments are required.\n\n"
    echo "$usage"
    exit 1
fi
exit 1
####CHECK BAM INPUT####
while read line
do
    if [ ! -e "$line" ]; then
         printf "Error. $line doesn't exist. Please check more carefully files passed in the 1st argument.\n"
         display_usage
         exit 1
    fi
done < "$1"
BAM_DATA=$(readlink -f "$1")

####SCRATCH WORKDIR####
OUTDIR=$(readlink -f "$2")
WORKDIR="/home/pedro.barbosa/scratch/statsBam"
if [ ! -d $WORKDIR ];then
    mkdir $WORKDIR
fi
if [ ! -d $OUTDIR ];then
    mkdir $OUTDIR
fi


cat > $WORKDIR/generateBamStats.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=getStatsFromBAM
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=10
#SBATCH --image=ummidock/bowtie2_samtools:latest
#SBATCH --workdir=$WORKDIR
#SBATCH --output=$WORKDIR/%j_getStats.log
#SBATCH --exclusive

srun="srun --exclusive -N1 -n1"
parallel="parallel --delay 0.2 -j \$SLURM_NTASKS --joblog parallel.log --resume"
time cat "$BAM_DATA" | \$parallel '\$srun shifter samtools view -bhF4 -@\$SLURM_CPUS_PER_TASK {} > $WORKDIR/{/.}.bam' 

#_stats.txt'

echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID 
#mv $WORKDIR/* $OUTDIR
EOL
sbatch $WORKDIR/generateBamStats.sbatch
