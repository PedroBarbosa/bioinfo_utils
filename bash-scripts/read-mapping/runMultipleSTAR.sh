#!/usr/bin/env bash
display_usage() {
echo 'Script to run STAR for each library. STAR executable must be in the path. If not, please edit STAR variable within the script with the full path.
Read groups are automatically added to each output based on the sample basename.
If there is annotation file available, please use it in the genome index generation step.
-1st argument must be the file listing RNA-seq pairs consecutively. One file per line.
-2nd argument must be the directory of the reference indexed database.
-3rd argument must be the number of threads to use.
-4th argument is optional. Flag to set 2-pass mappings. Recommended for alternative splicing analysis, more sensitive novel junction discovery. Available options: [true|false]. Default: "false".
-5th argument is optional. Do you want the output to be compatible with Cufflinks? If so, set "true". Available options: [true|false]. Default: false.
-6th argument is optional. Flag to keep alignments that contain non-canonical junctions. Available options: [true|false]. Default: "false", non-canonical junctions
are removed by default.
-7th argument is optional. Flag to output alginments in transcript coordinates compatible with eXpress software (allow indels and soft clippings in the
alignments). Available options: [true|false]. Default: "false".
'

}

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
        printf "Please provide at least the 3 first arguments required for the script.\n\n"
        display_usage
        exit 1
fi

STAR="STAR"
GENOMIC_INDEX="$2"
THREADS="$3"

###2-pass mappings###
if [ -z "$4" ] || [ "$4" = "false" ] ; then
    TWO_PASS_MAPPING="None"
elif [ "$4" = "true" ]; then
    TWO_PASS_MAPPING="Basic"
else
    printf "Please set a valid value for the 4th argument."
    display_usage
    exit 1
fi

##Cufflinks compatibility###
if [ -z "$5" ] || [ "$5" = "false" ] ; then
    CUFFLINKS_COMPATIBLE="None"
elif [ "$5" = "true" ]; then
    CUFFLINKS_COMPATIBLE="intronMotif"
else
    printf "Please set a valid value for the 5th argument."
    display_usage
    exit 1
fi

###Remove canonical junctions###
if [ -z "$6" ] || [ "$6" = "false" ] ; then
    NON_CANONICAL_JUNCTIONS="RemoveNoncanonical"
elif [ "$6" = "true" ]; then
    NON_CANONICAL_JUNCTIONS="None"
else
    printf "Please set a valid value for the 6th argument."
    display_usage
    exit 1
fi

###Express transcriptome aligments compatibility###
if [ -z "$7" ] || [ "$7" = "false" ] ; then
    EXPRESS_COMPATIBLE="RemoveNoncanonical"
elif [ "$7" = "true" ]; then
    EXPRESS_COMPATIBLE="Singleend"
else
    printf "Please set a valid value for the 7th argument."
    display_usage
    exit 1
fi

first_pair=true
second_pair=false


##SOME USEFUL PARAMETER SETTINGS THAT MIGHT PRESENT BETTER RESULTS:
#--outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --outFilterMatchNmin 60

while read line
do
        filename=$line
        bam_basename=$(basename $filename | sed 's/\(.*\)_.*/\1/')
        #SAMPLE_NAME=$(echo $pair1 | rev | cut -f 1 -d "/" | rev | cut -f1,2,3,4 -d "_")

        #add paired end pairs to the command
        if [ -f "$filename" -a "$first_pair" = "true" ]; then
                pair1=$filename
                first_pair=false
                second_pair=true

        elif [ -f "$filename" -a  "$second_pair" = "true" ]; then
                pair2=$filename
                first_pair="true"
                second_pair="false"
                printf "Running STAR for $bam_basename sample...\n"
                RUN_STAR="$STAR --runThreadN $THREADS --runMode alignReads --genomeDir $GENOMIC_INDEX --readFilesIn $pair1 $pair2 --outFileNamePrefix $bam_basename \
                  --outSAMattributes All --outSAMstrandField $CUFFLINKS_COMPATIBLE --outFilterIntronMotifs $NON_CANONICAL_JUNCTIONS --outSAMtype BAM SortedByCoordinate \
                  --chimSegmentMin 20 --outReadsUnmapped Fastx --quantMode TranscriptomeSAM GeneCounts --quantTranscriptomeBan $EXPRESS_COMPATIBLE \
                   --twopassMode $TWO_PASS_MAPPING --outSAMattrRGline ID:${bam_basename}_ID SM:$bam_basename PL:illumina LB:lib"

                printf "Command used for $bam_basename:\n$RUN_STAR\n\n"
                $RUN_STAR
                printf "Mapping completed for $bam_basename sample...\n\n\n"
        elif [ ! -f "$filename"  ]; then
                printf "$filename" is not a file
                exit 1
        fi
done <$1