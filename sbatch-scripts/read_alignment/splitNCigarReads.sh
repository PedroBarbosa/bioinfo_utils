#!/bin/bash
display_usage(){
 printf "Script to automatically run GATK4 splitNcigarReads utility on a set of bam files.\n
 Usage:.
    -1st argument must be the file containing the sam/bam files to mark duplicates. One file per line.
    -2nd argument must be the name of the ouptut directory. If '-' is given, output will be written in the parent directory of each bam file.
    -3rd argument is optional. Refers to the genome sequence in fasta format. If argument is not set or it gets '-' value, reference human genome hg38 will be used. Default: '/home/pedro.barbosa/mcfonseca/shared/genomes/human/hg38/GRCh38.primary.genome.fa'.\n"
}

if [ -z "$1" ] || [ -z "$2" ]; then
    printf "Error. Please set the required parameters of the script\n"
    display_usage
    exit 1
fi

####CHECK BAM/SAM INPUT####
while read line
do
    if [ ! -e "$line" ]; then
         printf "Error. $line doesn't exist. Please check more carefully files passed in the 1st argument.\n"
         display_usage
         exit 1
    fi
done < "$1"
BAM_DATA=$(readlink -f "$1")
JOBS=$(cat $BAM_DATA | wc -l)


####SCRATCH WORKDIR####
WORKDIR="/home/pedro.barbosa/scratch/gatk/splitNcigarString"
if [ ! -d $WORKDIR ];then
    mkdir $WORKDIR
fi

####OUTDIR####
if [ "$2" = "-" ]; then
    OUTDIR="{//}"
else
    OUTDIR=$(readlink -f "$2")
    if [ ! -d $OUTDIR ];then
        mkdir $OUTDIR
    fi
fi

####FASTA REFERENCE####
if [[ -z "$3" || "$3" == "-" ]];then
    REF_GENOME="/home/pedro.barbosa/mcfonseca/shared/genomes/human/hg38/GRCh38.primary.genome.fa"
else
    REF_GENOME=$(readlink -f "$3")
fi

JAVA_Xmx="--java-options '-Xmx45G'"
CMD="gatk ${JAVA_Xmx} SplitNCigarReads --tmp-dir $WORKDIR --max-mismatches-in-overhang 2 -R $REF_GENOME"

cat > runSplitNCigarReads.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=gatk_splitNcigar
#SBATCH --array=0-$(( $JOBS - 1 ))%10
#SBATCH --time=72:00:00
#SBATCH --mem=120G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --image=docker:broadinstitute/gatk:latest
##SBATCH --workdir=$WORKDIR
#SBATCH --output=%j_gatk_splitNcigar.log

readarray -t bams < <(cat $BAM_DATA)

SCRATCH_OUTDIR="$WORKDIR/\$SLURM_JOB_ID"
mkdir \$SCRATCH_OUTDIR
cd \$SCRATCH_OUTDIR

OUT_BASENAME=\$(basename \${bams[\$SLURM_ARRAY_TASK_ID]} .bam)
OUT_BASENAME=\${OUT_BASENAME/Aligned.sortedByCoord.out.bam/}

srun shifter $CMD -I \${bams[\$SLURM_ARRAY_TASK_ID]} -O \${OUT_BASENAME}_splitN.bam

if [ $OUTDIR = "{//}" ]; then
    mv * \$(dirname \${bams[\$SLURM_ARRAY_TASK_ID]})
else
    mv * $OUTDIR
fi

cd ../ && rm -rf \$SLURM_JOB_ID
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID 
if [ $OUTDIR != "{//}" ]; then
    mv \$SCRATCH_OUTDIR/* $OUTDIR
    mv ..\$SCRATCH_OUTDIR*log $OUTDIR

else
    echo "Since you set the output directory to be different for each sample, the run logs weren't moved anywhere. You can find them in the $WORKDIR folder."
fi
echo "\$(date): All done."
EOL
sbatch runSplitNCigarReads.sbatch
