#!/usr/bin/env bash
display_usage() {
echo 'Script to run STAR for multiple fastq files.
Read groups are automatically added to each output based on the sample basename.

-1st argument must be the file listing RNA-seq pairs consecutively. One file per line.
-2nd argument must be the directory of the reference indexed database. If "-", default hg38 star index will be used. 
-3rd argument must be the output directory.
-4th argument is optional. It is the identifier to extract the sample pair names from fastq files. If "-" is set, defaults are employed. Default: "_1.fq".
-5th argument is optional. Flag to output file compatible with Cufflinks and StringTie. Please set this flag to true if your data in unstranded. Available options: [true|false]. Default: false.
-6th argument is optional. Flag to add additional set of parameters to STAR command regarding the cases when the RNA insert size is smaller than the read length*2. [true|false]. Default: false.
-7th argument is optional. Flag to trigger read pair merging if pairs overlap > 10bp. [true|false|-]. Default: false.
-8th argument is optional. Flag to output wiggle format. Available options: [true|false]. Default: false.
-9th argument is optional. Flag to set 2-pass mappings. Only usefull when 4th argument is true. Recommended for alternative splicing analysis, more sensitive novel junction discovery. Available options: [true|false]. Default: "false".
-10th argument is optional. Flag to keep alignments that contain non-canonical junctions. Available options: [true|false]. Default: "false", non-canonical junctions
are removed by default.
-11th argument is optional. Flag to output alginments in transcript coordinates compatible with eXpress software (allow indels and soft clippings in the
alignments). Available options: [true|false]. Default: "false".'
}


if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] ; then
        printf "ERROR:Please provide at least the 3 first arguments required for the script.\n\n"
        display_usage 
        exit 1
fi

#Fastq
FASTQ=$(readlink -f "$1")
readarray fastq_array < $(readlink -f "$1")
JOBS=$(( ${#fastq_array[@]} / 2 ))

##INDEX##
if [[ $2 == "-" ]]; then
    INDEX="/home/pedro.barbosa/mcfonseca/shared/genomes/human/hg38/star/"
else
    INDEX=$(readlink -f "$2")
fi

#OUTPUT DIR##
if [[ ! -d $(readlink -f "$3") ]]; then
    mkdir $(readlink -f "$3")
fi
OUT="$(readlink -f "$3")"

##PAIR STRING##
if [[ -z "$4" || "$4" == "-" ]]; then
    DEL="_1.fq"
    p2="_2.fq"
else
    DEL="$4"
    p2=${DEL/1/2}
fi


##Cufflinks compatibility###
if [ -z "$5" ] || [ "$5" = "false" ] ; then
    CUFFLINKS_COMPATIBLE="None"
elif [ "$5" = "true" ]; then
    CUFFLINKS_COMPATIBLE="intronMotif"
else
    printf "Please set a valid value for the 5th argument.\n"
    display_usage
    exit 1
fi

###Small insert###
if [ "$6" = "true" ]; then
    #--alignEndsProtrude 30 ConcordantPair
    SMALL_INSERT="--outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --outFilterMatchNmin 60"
elif [ -z "$6" ] || [ "$6" = "false" ]; then
    SMALL_INSERT=""
else
    printf "Please set a valid value for the 6th argument.\n"
    display_usage
    exit 1
fi

##pairs merging###
if [ "$7" == "true" ]; then
    PAIRS_MERGE="--peOverlapNbasesMin 10 --chimOutType WithinBAM"
elif [ -z "$7" ] || [ "$7" == "false" ] || [ "$7" == "-" ]; then
    PAIRS_MERGE="--peOverlapNbasesMin 0"
else
    printf "Please set a valid value for the 7th argument.\n"
    display_usage
    exit 1
fi

###wiggle###
if [ "$8" = "true" ]; then
    WIGGLE="--outWigType wiggle"
elif [ -z "$8" ] || [ "$8" = "false" ]; then
    WIGGLE="None"
else
    printf "Please set a valid value for the 8th argument.\n"
    display_usage
    exit 1
fi

###2-pass mappings###
if [ -z "$9" ] || [ "$9" = "false" ] ; then
    TWO_PASS_MAPPING="None"
elif [ "$9" = "true" ]; then
    TWO_PASS_MAPPING="Basic"
else
    printf "Please set a valid value for the 9th argument.\n"
    display_usage
    exit 1
fi

###Remove canonical junctions###
if [ -z "${10}" ] || [ "${10}" = "false" ] ; then
    NON_CANONICAL_JUNCTIONS="RemoveNoncanonical"
elif [ "${10}" = "true" ]; then
    NON_CANONICAL_JUNCTIONS="None"
else
    printf "Please set a valid value for the 10th argument.\n"
    display_usage
    exit 1
fi

###Express transcriptome aligments compatibility###
if [ -z "${11}" ] || [ "${11}" = "false" ] ; then
    EXPRESS_COMPATIBLE="IndelSoftclipSingleend"
elif [ "${11}" = "true" ]; then
    EXPRESS_COMPATIBLE="Singleend"
else
    printf "Please set a valid value for the 11th argument.\n"
    display_usage
    exit 1
fi

CMD="STAR --runThreadN 20 --readFilesCommand gunzip -c --runMode alignReads --genomeDir $INDEX --outSAMattributes All --outSAMstrandField $CUFFLINKS_COMPATIBLE --outFilterIntronMotifs $NON_CANONICAL_JUNCTIONS $WIGGLE --outSAMtype BAM SortedByCoordinate --chimSegmentMin 20 --outReadsUnmapped Fastx --quantMode GeneCounts --twopassMode $TWO_PASS_MAPPING $PAIRS_MERGE $SMALL_INSERT"

cat > runSTAR.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=rna_star
#SBATCH --time=72:00:00
#SBATCH --mem=240G
#SBATCH --array=0-$(($JOBS - 1))%5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --image=dceoy/star:latest
#SBATCH --output=%j_star.log

SCRATCH_OUTDIR="/home/pedro.barbosa/scratch/star/\$SLURM_JOB_ID"
mkdir \$SCRATCH_OUTDIR && cd \$SCRATCH_OUTDIR 
echo "`date`: Analysis started."
srun="srun  -N1 -n1"

readarray -t pair1 < <(cat \$(readlink -f "$FASTQ") | grep -e "R1.fastq" -e "R1.fq" -e "_1.fastq" -e "_1.fq")
dir=\$(dirname \${pair1[\$SLURM_ARRAY_TASK_ID]})
pair2=\$(basename \${pair1[\$SLURM_ARRAY_TASK_ID]/$DEL/$p2})
fullpathpair2="\$dir/\$pair2"

outbasename=\${fullpathpair2##*/}
outbasename=\${outbasename%.*}
outbasename=\${outbasename/$p2/}

ulimit -n 16384
\$srun shifter $CMD --readFilesIn \${pair1[\$SLURM_ARRAY_TASK_ID]} \$fullpathpair2 --outFileNamePrefix \$outbasename --outSAMattrRGline ID:\${outbasename}_id SM:\${outbasename} PL:illumina LB:lib
mv * $OUT

#parallel="parallel --tmpdir \$SCRATCH_OUTDIR --halt soon,fail=1 --delay 0.2 -j $NTASKS --joblog parallel.log --resume-failed"
#cat "$FASTQ" | grep -v "R2" | \$parallel 'ulimit -n 16384; \$srun shifter $CMD --readFilesIn {} {=s/R1/R2/=} --outFileNamePrefix {=s{.*/}{};s/\_[^_]+$//;=}_ --outSAMattrRGline ID:{=s{.*/}{};s/\_[^_]+$//;=}_id SM:{=s{.*/}{};s/\_[^_]+$//;=} PL:illumina LB:lib; mv {=s{.*/}{};s/\_[^_]+$//;=}* $OUT'
#single end
#cat "$FASTQ" | grep -v "R2" | \$parallel 'ulimit -n 16384; \$srun shifter $CMD --readFilesIn {} --outFileNamePrefix {=s{.*/}{};s/\_[^_]+$//;=}_ --outSAMattrRGline ID:{=s{.*/}{};s/\_[^_]+$//;=}_id SM:{=s{.*/}{};s/\_[^_]+$//;=} PL:illumina LB:lib; mv {=s{.*/}{};s/\_[^_]+$//;=}* $OUT'
#mv parallel* $OUT
#cd ../ && rm -rf \$SLURM_JOB_ID
EOL
sbatch runSTAR.sbatch
