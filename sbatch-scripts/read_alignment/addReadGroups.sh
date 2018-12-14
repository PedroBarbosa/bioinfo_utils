#!/bin/bash
display_usage() {
printf 'Script to run GATK4 utility in lobo to add or replace read group information in a set of bam file. It will also create a new bam index for the new generated files.\n
    1st argument - Must be the list of sam/bam files to process.
    2nd argument - Must be the output directory to write the files. Values: [-|DIR]. "-" will write the new files in the dirname of each bam file.
    3rd argument is optional. To skip this argument and use defaults, set "-" here. Refers to the required addOrReplaceReadGroups picard utility LB,Pl,PU values. SM and ID will be automatically extracted from the filename, which is the more important stuff for downstream analysis. By default, LB=pe,PL=hiseq,PU=4000. You can set these values by splitting each value with a ",". E.g: "LB=mp,PL=miseq,PU=latest
    4th argument is optional. It refers to the sort order of the output file. If "-" is set, output is written as the input came in. Values:[unsorted|queryname|coordinate|duplicate]\n'    
}
if [ -z "$1" ] || [ -z "$2" ]; then
    printf "ERROR. Please set the required arguments for the script.\n"
    display_usage
    exit 1
fi

##INPUT##
if [ -e "$1" ]; then
    while read line
    do
        if [ ! -e "$line" ]; then
            printf "ERROR. Bam file "$line" does not exist. Please check the list of sam/bam files written in the 1st argument.\n"
            display_usage
            exit 1
        elif [ ${line: -4} != ".bam" -a ${line: -4} != ".sam" ] ; then
            printf "ERROR. Files listed in the 1st argument must be in bam or sam format.\n"
            display_usage
            exit 1
        fi
    done < "$1"
else
    printf "ERROR. File listing a set of bams does not exist. Please set a valid file for the 1st argument.\n"
    display_usage
    exit 1 
fi
BAM_DATA=$(readlink -f "$1")

##OUTDIR##
if [ -d "$2" ]; then
    printf "Output directory does not exist. Please provide a valid one, or set the 2nd argument to "-".\n"
    display_usage
    exit 1
elif [ "$2" = "-" ]; then
    OUTDIR="{//}"
else
    OUTDIR=$(readlink -f "$2")
fi
WORKDIR="/home/pedro.barbosa/scratch/gatk"
if [ ! -d "$WORKDIR" ]; then
    mkdir "$WORKDIR"
fi
##READGROUP INFO##
READGROUP=""
if [ -z "$3" ] || [ "$3" = "-" ]; then
    READGROUP="-LB pe -PL hiseq -PU 4000"
else
    IFS=',' 
    read -r -a array <<< "$3"
    if [ ${#array[@]} = 3 ]; then    
        for elem in "${array[@]}"
        do
            if [[ $elem == *"LB="* || $elem == *"PL="* || $elem == *"PU="* ]]; then
                READGROUP="$READGROUP-${elem/=/ } "
            else
                printf "ERROR.Please set the proper sintax for the 3rd argument. Only LB,PL and PU fields are allowed, and each of them must be assigned by a "=" character.\n"
                display_usage
                exit 1
            fi
        done
    else
        printf "ERROR. 3 fields are required for the 3rd argument (read groups). You set a different number.\n"
        display_usage
        exit 1
    fi
fi

##SORT ORDER OF THE OUTPUT##
SORT_ORDER=""
sort_values=("unsorted" "queryname" "coordinate" "duplicate")
if [ -z "$4" ] || [ "$4" == "-" ]; then
    printf "INFO. Sort order of the ouptut files will be kept.\n"
elif [[ !  " ${sort_values[@]} " =~ " ${4} " ]]; then
    printf "Please set a valid value for the 4th argument.\n"
    display_usage
    exit 1
else
    SORT_ORDER="-SO ${4} "
fi
CMD="gatk AddOrReplaceReadGroups ${SORT_ORDER}--CREATE_INDEX true $READGROUP"
cat > $WORKDIR/addReadGroups.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=gatk_rg
#SBATCH --nodes=2
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=2
#SBATCH --mem=230G
#SBATCH --time=72:00:00
#SBATCH --output=$WORKDIR/%j_gatk_rg.log
#SBATCH --image=docker:broadinstitute/gatk:latest
#SBATCH --workdir=$WORKDIR

SCRATCH_OUTDIR=$WORKDIR/\$SLURM_JOB_ID
if [ ! -d \$SCRATCH_OUTDIR ]; then
    mkdir \$SCRATCH_OUTDIR
fi
cd \$SCRATCH_OUTDIR

timestamp() {
    date +"%Y-%m-%d  %T"
}
export -f timestamp
srun="srun -N1 -n1 -c2"
parallel="parallel --delay 0.2 -j \$SLURM_NTASKS --env timestamp --joblog parallel.log --resume-failed"
echo "\$(timestamp) -> Analysis started!"
cat $BAM_DATA | \$parallel '\$srun shifter $CMD -SM {=s{.*/}{};s/\_[^_]+$//;s/\_[^_]+$//;s/\_[^_]+$//;=} -ID {=s{.*/}{};s/\_[^_]+$//;s/\_[^_]+$//;s/\_[^_]+$//;=}_id -I {} -O {/.}_RG.bam ; mv {/.}_RG.bam* $OUTDIR; echo -e "\$(timestamp) -> Finished parallel job number {#} (sample {/.})"'
echo -e "\$(timestamp) -> All done!"
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
EOL
sbatch $WORKDIR/addReadGroups.sbatch
sleep 2
cd $WORKDIR
mv addReadGroups.sbatch $(ls -td -- */ | head -n 1) 
