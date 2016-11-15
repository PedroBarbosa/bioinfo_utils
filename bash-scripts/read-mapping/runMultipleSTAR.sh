#!/usr/bin/env bash

display_usage() {
echo 'Script to run STAR for each library. STAR executable must be in the path. If not, please edit STAR variable within the script with the full path.
Read groups are automatically added to each output based on the sample basename.
If there is annotation file available, please use it in the genome index generation step [check 4th argument].


-1st argument must be the file listing RNA-seq pairs consecutively. One file per line.
-2nd argument must be the directory of the reference indexed database.
-3rd argument must be the number of threads to use.
-4th argument must be flag indicating if there is annotation available. Available option: [true|false]. Default: yes,expected to be added in the index generation step. Argument is required to know
if parameteres regarding gene/transcripts quantification should be on the STAR command. If no annotation available, 9th argument is useless for instance. Please check STAR manual for further details.
-5th argument is optional. Flag to add additional set of parameters to STAR command regarding the cases when the RNA insert size is smaller than the read length*2. [true|false]. Default: false.
-6th argument is optional. Flag to output wiggle format. Available options: [true|false]. Default: false.
-7th argument is optional. Flag to set 2-pass mappings. Only usefull when 4th argument is true. Recommended for alternative splicing analysis, more sensitive novel junction discovery. Available options: [true|false]. Default: "false".
-8th argument is optional. Flag to output file compatible with Cufflinks and StringTie. Please set this flag to true if your data in unstranded. Available options: [true|false]. Default: false.
-9th argument is optional. Flag to keep alignments that contain non-canonical junctions. Available options: [true|false]. Default: "false", non-canonical junctions
are removed by default.
-10th argument is optional. Flag to output alginments in transcript coordinates compatible with eXpress software (allow indels and soft clippings in the
alignments). Available options: [true|false]. Default: "false".'
}


if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] ; then
        printf "ERROR:Please provide at least the 4 first arguments required for the script.\n\n"
        display_usage
        exit 1
fi

STAR="STAR"
GENOMIC_INDEX="$2"
THREADS="$3"

###annotation available###
if [ "$4" != "true" ] && [ "$4" != "false" ] ; then
    printf "Please set a valid value for the 4th argument.\n"
    display_usage
    exit 1
fi

###Small insert###
if [ "$5" = "true" ]; then
    SMALL_INSERT="--outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --outFilterMatchNmin 60"
elif [ -z "$6" ] || [ "$6" = "false" ]; then
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

##Cufflinks compatibility###
if [ -z "$8" ] || [ "$8" = "false" ] ; then
    CUFFLINKS_COMPATIBLE="None"
elif [ "$8" = "true" ]; then
    CUFFLINKS_COMPATIBLE="intronMotif"
else
    printf "Please set a valid value for the 8th argument.\n"
    display_usage
    exit 1
fi

###Remove canonical junctions###
if [ -z "$9" ] || [ "$9" = "false" ] ; then
    NON_CANONICAL_JUNCTIONS="RemoveNoncanonical"
elif [ "$9" = "true" ]; then
    NON_CANONICAL_JUNCTIONS="None"
else
    printf "Please set a valid value for the 9th argument.\n"
    display_usage
    exit 1
fi

###Express transcriptome aligments compatibility###
if [ -z "${10}" ] || [ "${10}" = "false" ] ; then
    EXPRESS_COMPATIBLE="IndelSoftclipSingleend"
elif [ "${10}" = "true" ]; then
    EXPRESS_COMPATIBLE="Singleend"
else
    printf "Please set a valid value for the 10th argument.\n"
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

		if [ "$4" = "true" ]; then
                	RUN_STAR="$STAR --runThreadN $THREADS --runMode alignReads --genomeDir $GENOMIC_INDEX --readFilesIn $pair1 $pair2 --outFileNamePrefix $bam_basename
 --outSAMattributes All --outSAMstrandField $CUFFLINKS_COMPATIBLE --outFilterIntronMotifs $NON_CANONICAL_JUNCTIONS $WIGGLE --outSAMtype BAM SortedByCoordinate
 --chimSegmentMin 20 --outReadsUnmapped Fastx --quantMode TranscriptomeSAM GeneCounts --quantTranscriptomeBan $EXPRESS_COMPATIBLE --twopassMode $TWO_PASS_MAPPING
 --outSAMattrRGline ID:${bam_basename}_ID SM:$bam_basename PL:illumina LB:lib $SMALL_INSERT"
		else
			RUN_STAR="$STAR --runThreadN $THREADS --runMode alignReads --genomeDir $GENOMIC_INDEX --readFilesIn $pair1 $pair2 --outFileNamePrefix $bam_basename
 --outSAMattributes All --outSAMstrandField $CUFFLINKS_COMPATIBLE --outFilterIntronMotifs $NON_CANONICAL_JUNCTIONS $WIGGLE --outSAMtype BAM SortedByCoordinate
 --chimSegmentMin 20 --outReadsUnmapped Fastx --twopassMode $TWO_PASS_MAPPING --outSAMattrRGline ID:${bam_basename}_ID SM:$bam_basename PL:illumina LB:lib $SMALL_INSERT"
		fi
                printf "Command used for $bam_basename:\n$RUN_STAR\n\n"
                $RUN_STAR
                printf "Mapping completed for $bam_basename sample...\n\n\n"
        elif [ ! -f "$filename"  ]; then
                printf "$filename" is not a file
                exit 1
        fi
done <$1
