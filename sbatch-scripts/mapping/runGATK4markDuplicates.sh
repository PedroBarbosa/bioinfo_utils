#!/bin/bash
display_usage(){
 printf "Script to automatically run GATK4 mark duplicates utility on a set of bam files.\n
 Usage:.
    -1st argument must be the file containing the sam/bam files to mark duplicates. One file per line.
    -2nd argument must be the name of the ouptut directory. If '-' is given, output will be written in the parent directory of each bam file.
    -3rd argument is optional. If set to true, the tool will remove duplicates instead of marking them in the output file. If '-' is set, this parameter will be ignored, and defaults will be employed (false).
    -4th argument is optional. If set to false, GNU parallel will be disabled to run the set of samples provided in the 1st argument. Options: [true|false]. Default: true, GNU parallel is used to parallelize the job.\n"
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

####SCRATCH WORKDIR####
WORKDIR="/home/pedro.barbosa/scratch/gatk"
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

CMD="/gatk/gatk-launch MarkDuplicatesWithMateCigar --REMOVE_DUPLICATES"
####INNDEX####
if [ -z "$3" ] || [ "$3" = "-" ] || [ "$3" = "false" ]; then
    CMD="$CMD false"
elif [ "$3" = "true" ]; then
    CMD="$CMD true"
else
    printf "Error. Please set a valid value for the third argument. [true|false|-]\n"
    display_usage
    exit 1    
fi

###PARALLEL####
if [ -z "$4" ] || [ "$4" = "true" ]; then
    NODES=1 #2 #1
    NTASKS=10 #10 #4
    CPUS=4 #8 #10
    PARALLEL=true
elif [ "$4" = "false" ]; then
    NODES=1
    NTASKS=1
    CPUS=40
    PARALLEL=false
else
    printf "Please set a valid value for the 4th argument\n"
    display_usage
    exit 1
fi

cat > $WORKDIR/runGATKmarkDup.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=gatk_md
#SBATCH --time=72:00:00
#SBATCH --mem=150G
#SBATCH --nodes=$NODES
#SBATCH --ntasks=$NTASKS
#SBATCH --cpus-per-task=$CPUS
#SBATCH --image=docker:broadinstitute/gatk:latest
#SBATCH --workdir=$WORKDIR
#SBATCH --output=$WORKDIR/%j_markDup.log

SCRATCH_OUTDIR="$WORKDIR/\$SLURM_JOB_ID"
mkdir \$SCRATCH_OUTDIR
cd \$SCRATCH_OUTDIR
if [ "$PARALLEL" = "true" ]; then
	srun="srun --exclusive -N1 -n1"
	parallel="parallel --tmpdir \$SCRATCH_OUTDIR --delay 0.2 -j $NTASKS --joblog parallel.log --resume"	
        cat "$BAM_DATA" | \$parallel '\$srun shifter $CMD -I={} -M={/.}_dup.metrics -O={/.}_DUP.bam; mv {/.}* $OUTDIR'
else
	for j in \$(find $BAM_DATA -exec cat {} \; );do
	    printf "\$(date): Processing \$(basename \$j) file.\n"
	    i=\$(basename \$j | rev | cut -d'.' -f2- | rev)
	    outbam="\${i}_DUP.bam
            metrics="\${i}_dup.metrics
	    echo "srun shifter $CMD -I=\$j -M=\$metrics -O=\$outbam"
            if [ $OUTDIR = "{//}" ]; then
                echo "mv \${i}* \$(dirname \$j)"
            else
                echo "mv \${i}* $OUTDIR"
            fi
            printf "\$(date): Done.\n"
	done
fi 
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID 
if [ $OUTDIR != "{//}" ]; then
    echo "mv \$SCRATCH_OUTDIR/* $OUTDIR"
else
    echo "Since you set the output directory to be different for each sample, the run logs weren't moved anywhere. You can find them in the $WORKDIR folder."
fi
echo "\$(date): All done."
EOL
sbatch $WORKDIR/runGATKmarkDup.sbatch



