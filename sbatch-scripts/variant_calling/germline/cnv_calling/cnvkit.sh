#!/bin/bash
display_usage(){
 printf "Script to automatically run cnvkit in lobo for multiple bam files.\n
 Usage:.
    -1st argument must be the file listing the tumor samples (1 bam per line).
    -2nd argument must be the file listing the normal samples (1 bam per line).
    -3rd argument must be the output directory.\n"
}

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
    printf "Error. Please set the required parameters of the script\n"
    display_usage
    exit 1
fi

####CHECK BAM INPUT####
tumor=""
normal="--normal"
while read line
do
    if [ ! -e "$line" ]; then
         printf "Error. $line doesn't exist. Please check more carefully files passed in the 1st argument.\n"
         display_usage
         exit 1
    else
        tumor="$tumor $line" 
    fi
done < "$(readlink -f "$1")"


while read line
do
    if [ ! -e "$line" ]; then
         printf "Error. $line doesn't exist. Please check more carefully files passed in the 2nd argument.\n"
         display_usage
         exit 1
    else
        normal="$normal $line"
    fi
done < "$(readlink -f "$2")"
    

####SCRATCH WORKDIR####
OUTDIR=$(readlink -f "$3")
WORKDIR="/home/pedro.barbosa/scratch/cnvkit"
if [ ! -d $WORKDIR ];then
    mkdir $WORKDIR
fi
if [ ! -d $OUTDIR ];then
    mkdir $OUTDIR
fi

####REFERENCE GENOME####
genome="/mnt/nfs/lobo/IMM-NFS/genomes/hg19/Sequence/WholeGenomeFasta/genome.fa"

cat > $WORKDIR/runCNVkit.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=cnvkit
#SBATCH --time=24:00:00
#SBATCH --mem=150G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --image=docker:etal/cnvkit
#SBATCH --workdir=$WORKDIR
#SBATCH --output=$WORKDIR/%j_cnvkit.log

SRUN="srun -N1 -n1 shifter"
cmd="cnvkit.py batch -m hybrid --count-reads --drop-low-coverage --scatter --diagram --fasta $genome"
\$SRUN \$cmd $normal -d $WORKDIR $tumor

echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID 
#mv $WORKDIR/* $OUTDIR
EOL
sbatch $WORKDIR/runCNVkit.sbatch
