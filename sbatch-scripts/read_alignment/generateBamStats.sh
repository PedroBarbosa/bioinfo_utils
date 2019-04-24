#!/bin/bash
display_usage(){
 printf "Script to automatically count alignments in a list of BAM files.\n
 Usage:.
    -1st argument must be the file containing the files sequences to count.
    -2nd argument must be the name of the ouptut directory.
    -3rd argument must be the name of the output file.
    -4th argument is optional.  Refers to the number of nodes,tasks and cpus per task, respectively, to employ on this slurm job in lobo (tasks will be set in parallel,not in the srun command). Default:2,8,10. '-' skips this argument.
    -5rd argument is optional. You may write here additional flags to pass on to samtools. Multiple flags may by setting the ';' delimiter.\n" 
}

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
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

##JOB SETTINGS""
if [[ -z "$4" || "$4" = "-" ]]; then
    NODES=2
    NTASKS=8
    CPUS_PER_TASK=10
else
    IFS=','
    read -r -a array <<< "$4"
    if [ ${#array[@]} = 3 ]; then
        for elem in "${array[@]}"
        do
            if ! [[ "$elem" =~ $re ]]; then
                printf "Error. Please set INT numbers for the number of nodes, tasks and cpus per task.\n"
                display_usage
                exit 1
            fi
        done
        NODES=${array[0]}
        NTASKS=${array[1]}
        CPUS_PER_TASK=${array[2]}
    else
        printf "ERROR. 3 fields are required for the 4th argument (nodes,tasks,cpus per task). You set a different number.\n"
        display_usage
        exit 1
    fi
fi

OUT_FILE="$3"
MAPPER="STAR"
cat > $WORKDIR/generateBamStats.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=getStatsFromBAM
#SBATCH --time=72:00:00
#SBATCH --mem=245G
#SBATCH --nodes=$NODES
#SBATCH --ntasks=$NTASKS
#SBATCH --cpus-per-task=$CPUS_PER_TASK
#SBATCH --image=ummidock/bowtie2_samtools:latest
##SBATCH --workdir=$WORKDIR
#SBATCH --output=%j_getStats.log

timestamp() {
    date +"%Y-%m-%d  %T"
}
#metrics="{=s{.*/}{};s/\_[^_]+$//;s/\_[^_]+$//;s/\_[^_]+$//=}\t\$reads\t\$aln\t\$unmapped\t\$duplicates\t\$primary_q10\t\$proper\t\$secondary\t\$chimeric\t\$unique"; \
export -f timestamp
scratch_out=$WORKDIR/\$SLURM_JOB_ID
mkdir \$scratch_out
cd \$scratch_out
srun="srun -N1 -n1 --slurmd-debug 4"
parallel="parallel -k --delay 0.2 -j \$SLURM_NTASKS  --env timestamp --joblog parallel.log --resume-failed"
echo "\$(timestamp) -> Analysis started!"
header="sample\t#reads\t#alignments\t#unmapped_reads\t#duplicate_reads\t#primary_linear_q10\t#proper_pair\t#secondary\t#chimeric\t#unique"
echo -e "\$header" > $OUT_FILE
time cat "$BAM_DATA" | \$parallel 'case "${MAPPER}" in "bwa-mem") unique=\$(\$srun shifter samtools view -F2308q10 -@\$SLURM_CPUS_PER_TASK {} | grep -cv "XA:") ;; "bowtie2") unique=\$(\$srun shifter samtools view -F2308q10 -@\$SLURM_CPUS_PER_TASK {} | grep -cv "XS:") ;; "STAR") unique=\$(\$srun shifter samtools view -cq255 -@\$SLURM_CPUS_PER_TASK {}) ;; "hisat2") unique=\$(\$srun shifter samtools view -@\$SLURM_CPUS_PER_TASK {} | grep -c "NH:i:1") ;; esac && \
\$srun shifter samtools sort -n -O bam -o {/.}_tmp.bam -@\$SLURM_CPUS_PER_TASK {} && \
fragments=\$(\$srun shifter samtools view {/.}_tmp.bam | cut -f1 | uniq | wc -l) && reads=\$((\$fragments * 2)) &&  rm {/.}_tmp.bam && \
aln=\$(\$srun shifter samtools view -cF4 -@\$SLURM_CPUS_PER_TASK {}) && \
unmapped=\$(\$srun shifter samtools view -cf4 -@\$SLURM_CPUS_PER_TASK {}) && \
duplicates=\$(\$srun shifter samtools view -cf1024 -@\$SLURM_CPUS_PER_TASK {}) && \
primary_q10=\$(\$srun shifter samtools view -cF2308q10 -@\$SLURM_CPUS_PER_TASK {}) && \
proper=\$(\$srun shifter samtools view -cf2 -@\$SLURM_CPUS_PER_TASK {}) && \
secondary=\$(\$srun shifter samtools view -cf256 -@\$SLURM_CPUS_PER_TASK {}) && \
chimeric=\$(\$srun shifter samtools view -cf2048 -@\$SLURM_CPUS_PER_TASK {}) && \
metrics="{/.}\t\$reads\t\$aln\t\$unmapped\t\$duplicates\t\$primary_q10\t\$proper\t\$secondary\t\$chimeric\t\$unique"; \
echo -e "\$metrics" >> $OUT_FILE && \
echo -e "\$(timestamp) -> Finished parallel job number {#} (sample {=s{.*/}{};s/\_[^_]+$//;s/\_[^_]+$//=})"'
mv $OUT_FILE $OUTDIR
mv ../\$SLURM_JOB_ID_getStats.log $OUTDIR
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
echo -e "\$(timestamp) -> All done!"
EOL
sbatch $WORKDIR/generateBamStats.sbatch
sleep 3 
cd $WORKDIR
mv generateBamStats.sbatch $(ls -td -- */ | head -n 1) 
