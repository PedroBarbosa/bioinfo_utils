#!/bin/bash
display_usage(){
 printf "Script to automatically convert a list of bam files into fastq. It will use bamtools for the job.\n
 Usage:
    -1st argument must be the list of BAM files to process.
    -2nd argument is optional. Refers to the format to convert the output. [default: fastq. Options:bed|fasta|fastq|json|pileup|sam|yaml]
    -2nd argument is optional. Referes to the output directory to store all fastq files. [default: output is written on the bam directory.\n"
}

if [ -z "$1" ]; then
    printf "Error. Please set the required parameters\n"
    display_usage
    exit 1
fi

while read line
do
    if [ ! -e "$line" ]; then
        printf "Error. $line doesn't exist. Please check more carefully files passed in the 1st argument.\n"
        display_usage
        exit 1
    fi
done < "$1"
bamfiles="$1"
workdir="/home/pedro.barbosa/scratch/bamConvert"
if [ ! -d $workdir ]; then
    mkdir $workdir
fi

same_out_dir=false
if [ -d "$2" ];
    echo "$2 directory will be used to write the final files"
    same_out_dir=true
elif [ -d "$2" ];
    printf "$2 directory doesn't exist. Please set a valid one!\n"
    display_usage
    exit 1

cat > $workdir/bam2fastq.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=hcm_convert
#SBATCH --output=/home/pedro.barbosa/scratch/bamConvert/%j_convert.log
#SBATCH --mail-user=pedro.barbosa@medicina.ulisboa.pt
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --workdir=$workdir
#SBATCH --image=mcfonsecalab/htstools_plus

if [ $same_out_dir = "false" ]; then

    cat "$FASTQ_DATA" | grep -v "R2" | \$parallel '\$srun shifter bwa mem -t $CPUS $INDEX {} {=s/R1/R2/=} | shifter --image=ummidock/bowtie2_samtools:latest samtools sort -O bam -@ $CPUS -o $WORKDIR/{/.}_sort.bam -; shifter --image=ummidock/bowtie2_samtools:latest samtools index $WORKDIR/{/.}_sort.bam'
srun="srun --exclusive -N1 -n1"
parallel="parallel --delay 0.2 -j \$SLURM_NTASKS --resume" #--joblog parallel.log"
cat $bamfiles | \$parallel --dryrun '\$srun shifter bamtools convert 

echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID 
EOL

sbatch $workdir/bam2fastq.sbatch
