#!/bin/bash

display_usage(){
 printf "Script to automatically count alignments in a list of BAM files.\n
 Usage:.
    -1st argument must be the file containing the files sequences to count.
    -2nd argument must be the name of the ouptut directory.
    -3rd argument is optional. You may write here additional flags to pass on to samtools. Multiple flags may by setting the ';' delimiter.\n" 
}

if [ -z "$1" ] || [ -z "$2" ]; then
    printf "Error. Please set the required parameters of the script\n"
    display_usage
    exit 1
fi

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

OUT_FILE="mappingStats.tsv"
MAPPER="bwa-mem"
cat > $WORKDIR/generateBamStats.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=getStatsFromBAM
#SBATCH --time=72:00:00
#SBATCH --mem=120G
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=8
#SBATCH --image=ummidock/bowtie2_samtools:latest
#SBATCH --workdir=$WORKDIR
#SBATCH --output=$WORKDIR/%j_getStats.log

timestamp() {
    date +"%Y-%m-%d  %T"
}
export -f timestamp
scratch_out=$WORKDIR/\$SLURM_JOB_ID
mkdir \$scratch_out
cd \$scratch_out
srun="srun -N1 -n1"
parallel="parallel -k --delay 0.2 -j \$SLURM_NTASKS  --env timestamp --joblog parallel.log --resume-failed"
echo "\$(timestamp) -> Analysis started!"
header="sample\t#alignments\t#unmapped_reads\t#duplicate_reads\t#primary_linear_q10\t#proper_pair\t#secondary\t#chimeric\t#unique"
echo -e "\$header" > $OUT_FILE
time cat "$BAM_DATA" | \$parallel 'case "${MAPPER}" in "bwa-mem") unique=\$(\$srun shifter samtools view -F2308q10 -@\$SLURM_CPUS_PER_TASK {} | grep -cv "XA:") ;; "bowtie2") unique=\$(\$srun shifter samtools view -F2308q10 -@\$SLURM_CPUS_PER_TASK {} | grep -cv "XS:") ;; "STAR") unique=\$(\$srun shifter samtools view -cq255 -@\$SLURM_CPUS_PER_TASK {}) ;; "hisat2") unique=\$(\$srun shifter samtools view -@\$SLURM_CPUS_PER_TASK {} | grep -c "NH:i:1") ;; esac && \
aln=\$(\$srun shifter samtools view -cF4 -@\$SLURM_CPUS_PER_TASK {}) && \
unmapped=\$(\$srun shifter samtools view -cf4 -@\$SLURM_CPUS_PER_TASK {}) && \
duplicates=\$(\$srun shifter samtools view -cf1024 -@\$SLURM_CPUS_PER_TASK {}) && \
primary_q10=\$(\$srun shifter samtools view -cF2308q10 -@\$SLURM_CPUS_PER_TASK {}) && \
proper=\$(\$srun shifter samtools view -cf2 -@\$SLURM_CPUS_PER_TASK {}) && \
secondary=\$(\$srun shifter samtools view -cf256 -@\$SLURM_CPUS_PER_TASK {}) && \
chimeric=\$(\$srun shifter samtools view -cf2048 -@\$SLURM_CPUS_PER_TASK {}) && \
metrics="{=s{.*/}{};s/\_[^_]+$//;s/\_[^_]+$//=}\t\$aln\t\$unmapped\t\$duplicates\t\$primary_q10\t\$proper\t\$secondary\t\$chimeric\t\$unique"; \
echo -e "\$metrics" >> $OUT_FILE && \
echo -e "\$(timestamp) -> Finished parallel job number {#} (sample {=s{.*/}{};s/\_[^_]+$//;s/\_[^_]+$//=})"'
mv $OUT_FILE $OUTDIR
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
echo -e "\$(timestamp) -> All done!"
EOL
sbatch $WORKDIR/generateBamStats.sbatch
sleep 1 
cd $WORKDIR
mv generateBamStats.sbatch $(ls -td -- */ | head -n 1) 
