#!/usr/bin/env bash
display_usage() {
echo 'Script to run STAR for multiple fastq files.
Read groups are automatically added to each output based on the sample basename.

-1st argument must be the file listing RNA-seq pairs consecutively. One file per line.
-2nd argument must be the directory of the reference indexed database.
-3rd argument must be the number of nodes,tasks and threads to use with GNU parallel (split numbers by ",").
-4th argument is optional. Flag to output file compatible with Cufflinks and StringTie. Please set this flag to true if your data in unstranded. Available options: [true|false]. Default: false.
-5th argument is optional. Flag to add additional set of parameters to STAR command regarding the cases when the RNA insert size is smaller than the read length*2. [true|false]. Default: false.
-6th argument is optional. Flag to output wiggle format. Available options: [true|false]. Default: false.
-7th argument is optional. Flag to set 2-pass mappings. Only usefull when 4th argument is true. Recommended for alternative splicing analysis, more sensitive novel junction discovery. Available options: [true|false]. Default: "false".
-8th argument is optional. Flag to keep alignments that contain non-canonical junctions. Available options: [true|false]. Default: "false", non-canonical junctions
are removed by default.
-9th argument is optional. Flag to output alginments in transcript coordinates compatible with eXpress software (allow indels and soft clippings in the
alignments). Available options: [true|false]. Default: "false".'
}


if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] ; then
        printf "ERROR:Please provide at least the 3 first arguments required for the script.\n\n"
        display_usage 
        exit 1
fi

#Fastq
FASTQ=$(readlink -f "$1")
INDEX=$(readlink -f "$2")
OUT=$PWD

#Resources
re='^[0-9]+$'
IFS=','
read -r -a array <<< "$3"
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
    if [ "$5" = "false" ]; then
        PARALLEL="false"
    elif [ "$5" != "true" ]; then
        printf "Error. Please set a valid value for the 5th argument (gnu parallel flag).\n"
        display_usage
        exit 1
    fi
else
    printf "ERROR. 3 fields are required for the 9th argument (nodes,tasks,cpus per task). You set a different number.\n"
    display_usage
    exit 1
fi

##Cufflinks compatibility###
if [ -z "$4" ] || [ "$4" = "false" ] ; then
    CUFFLINKS_COMPATIBLE="None"
elif [ "$4" = "true" ]; then
    CUFFLINKS_COMPATIBLE="intronMotif"
else
    printf "Please set a valid value for the 4th argument.\n"
    display_usage
    exit 1
fi

###Small insert###
if [ "$5" = "true" ]; then
    SMALL_INSERT="--outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --outFilterMatchNmin 60"
elif [ -z "$5" ] || [ "$5" = "false" ]; then
    SMALL_INSERT=""
else
    printf "Please set a valid value for the 5th argument.\n"
    display_usage
    exit 1
fi

###wiggle###
if [ "$6" = "true" ]; then
    WIGGLE="--outWigType wiggle"
elif [ -z "$6" ] || [ "$6" = "false" ]; then
    WIGGLE="None"
else
    printf "Please set a valid value for the 6th argument.\n"
    display_usage
    exit 1
fi


###2-pass mappings###
if [ -z "$7" ] || [ "$7" = "false" ] ; then
    TWO_PASS_MAPPING="None"
elif [ "$7" = "true" ]; then
    TWO_PASS_MAPPING="Basic"
else
    printf "Please set a valid value for the 7th argument.\n"
    display_usage
    exit 1
fi

###Remove canonical junctions###
if [ -z "$8" ] || [ "$8" = "false" ] ; then
    NON_CANONICAL_JUNCTIONS="RemoveNoncanonical"
elif [ "$8" = "true" ]; then
    NON_CANONICAL_JUNCTIONS="None"
else
    printf "Please set a valid value for the 8th argument.\n"
    display_usage
    exit 1
fi

###Express transcriptome aligments compatibility###
if [ -z "$9" ] || [ "$9" = "false" ] ; then
    EXPRESS_COMPATIBLE="IndelSoftclipSingleend"
elif [ "$9" = "true" ]; then
    EXPRESS_COMPATIBLE="Singleend"
else
    printf "Please set a valid value for the 9th argument.\n"
    display_usage
    exit 1
fi

CMD="STAR --runThreadN $CPUS_PER_TASK --readFilesCommand gunzip -c --runMode alignReads --genomeDir $INDEX --outSAMattributes All --outSAMstrandField $CUFFLINKS_COMPATIBLE --outFilterIntronMotifs $NON_CANONICAL_JUNCTIONS $WIGGLE --outSAMtype BAM SortedByCoordinate --chimSegmentMin 20 --outReadsUnmapped Fastx --quantMode GeneCounts --twopassMode $TWO_PASS_MAPPING $SMALL_INSERT"

cat > runSTAR.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=rna_star
#SBATCH --time=72:00:00
#SBATCH --mem=240G
#SBATCH --nodes=$NODES
#SBATCH --ntasks=$NTASKS
#SBATCH --cpus-per-task=$CPUS_PER_TASK
#SBATCH --image=argrosso/star:2.6.0b
#SBATCH --output=%j_star.log

SCRATCH_OUTDIR="/home/pedro.barbosa/scratch/star/\$SLURM_JOB_ID"
mkdir \$SCRATCH_OUTDIR && cd \$SCRATCH_OUTDIR 
echo "`date`: Analysis started."
srun="srun --exclusive -N1 -n1"
parallel="parallel --tmpdir \$SCRATCH_OUTDIR --halt soon,fail=1 --delay 0.2 -j $NTASKS --joblog parallel.log --resume-failed"
cat "$FASTQ" | grep -v "R2" | \$parallel 'ulimit -n 2048; \$srun shifter $CMD --readFilesIn {} {=s/R1/R2/=} --outFileNamePrefix {=s{.*/}{};s/\_[^_]+$//;=}_ --outSAMattrRGline ID:{=s{.*/}{};s/\_[^_]+$//;=}_id SM:{=s{.*/}{};s/\_[^_]+$//;=} PL:illumina LB:lib; mv {=s{.*/}{};s/\_[^_]+$//;=}* $OUT'
cd ../ && rm -rf \$SLURM_JOB_ID
EOL
sbatch runSTAR.sbatch