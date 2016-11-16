#!/bin/bash

display_usage() {
echo 'Script to run multiple StringTie for several bam files. Stringtie must be in the system path. If not, please set full location within the script.

-1st argument must be the file listing the BAM files to perform transcriptome assembly. One file per line. If 8th argument is set to true, this file should represent a set of GTF/GFF files to merge.
-2nd argument must be the path for reference annotation (GTF or GFF) for guided transcriptome assembly. If there is no annotation, please set "-" on this argument. Available options: [FILE|-].
-3nd argument must be the number of threads to run StringTie.
-4th argument must be the output directory. If not exists, will create automatically.
-5rd argument is optional. Flag if one wants to run StringTie with stringent settings for transcript identification. Available options: [true|false]. Default:false. If true, parameteres "-c", "-j", "-a", "-m", "-f" will be affected.
-6th argument is optional. Output tables required for Ballgown downstream analysis. Requires reference annotation provided (2nd argument). Available options: [true|false]. Default:false.
-7th argument is optional. Output only assembled transcripts that match reference transcripts. Requires reference annotation provided (2nd argument.). Usefull for a 2nd round of StringTIe.
-8th argument is optional. Flag to activate TRANSCRIPT MERGE MODE. Available options: [true|false]. Default: false (normal transcriptome assembly). If set, it will merge a set of input GTF/GFF files provided in the 1st argument into a non redundant set of transcripts. If 2nd argument is set, merging will include the reference, otherwise (if "-" set), merge is done based only on the list of GTF/GFF provided.'

}

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] ; then
	printf "Error:Please provide the required arguments for the scripts.\n\n"
	display_usage
	exit 1
fi

STRINGTIE="stringtie"
THREADS="-p $3"
if [ -d "$4" ]; then
    OUTPUT="-o $4"
else
    mkdir "$4"
    OUTPUT="-o $4"
fi

CMD="$STRINGTIE $THREADS"
#####################
if [ -f "$2" ]; then
    REFERENCE_ANNOTATION="-G $2"
    CMD="$CMD $REFERENCE_ANNOTATION"
elif [ "$2" != "-" ]; then
    printf "Please provide a valid option for the 2nd argument.\n\n"
    display_usage
    exit 1
fi
#####################
#####################
if [ "$5" = "true" ]; then
    STRINGENT="-c 5 -j 1.5 -a 20 -m 400 -f 0.2"
    CMD="$CMD $STRINGENT"
elif [ "$5" != "false" ]; then
    printf "Please provide a valid option for the 5th argument.\n\n"
    display_usage
    exit 1
fi
######################
######################
if [ "$6" = "true" -a "$2" = "-" ]; then
    printf "Ballgown tables can only be outputed when reference annotation is provided (2nd argument). Please set false for this 6h argument or provide a valid file in the 2nd argument.\n\n"
    display_usage
    exit 1
elif [ -f "$2" -a "$6" = "true" ]; then
    CMD="$CMD -B"
elif [ "$6" != "false" ]; then
    printf "Please provide a valid option for the 6th argument.\n\n"
    display_usage
    exit 1
fi
######################
######################
if [ "$7" = "true" -a "$2" = "-" ]; then
    printf "7th argument is only valid when reference annotation is provided (2nd argument).  Please set false for this 7h argument or provide a valid file in the 2nd argument.\n\n"
    display_usage
    exit 1
elif [ -f "$2" -a "$7" = "true" ]; then
    CMD="$CMD -e"
elif [ "$7" != "false" ]; then
    printf "Please provide a valid option for the 7th argument.\n\n"
    display_usage
    exit 1
fi
#####################
#####################
if [ "$8" = "true" -a "$2" = "-" ]; then
    CMD="$STRINGTIE --merge $THREADS "
elif [ -f "$2" -a "$8" = "true" ]; then
    CMD="$STRINGTIE --merge $THREADS $REFERENCE_ANNOTATION"
elif [ "$8" != "false" ]; then
    printf "Please provide a valid option for the 8th argument.\n\n"
    display_usage
    exit 1
fi
#####################


########RUNNING################
if [ "$8" = "true" ]; then
    printf "Running StringTie to merge individual GTF files..\n"
    while read line
    do
        if [[ $line =~ \.gtf$ ]] ; then
            continue
        else
            printf "'.gtf' extension not found in at least one of the input files. Please check if you are providing a list of bam files instead of gtf. gtf extension is required to run in the merge mode.\n"
            display_usage
            exit 1
        fi
    done < $1

    $CMD="$CMD $OUTPUT/stringTie_merged.gtf $1"
    printf "$CMD"

else
    printf "Running StringTie to assemble transcript from BAM files..\n"
    while read line
    do
	    FILENAME=$line
	    #BASENAME=$(basename ${FILENAME} .bam)

	    #STAR SPECIFIC
	    BASENAME=$(basename ${FILENAME} .sortedByCoord.out.bam)
	    printf "Started sample $BASENAME...\n"

	    #Add coverage tables if annotation provided
	    if [ -f "$2" ]; then
	        $CMD="$CMD -o $OUTPUT/$BASENAME-transcripts.gtf -A $OUTPUT/$BASENAME-geneAbundances.txt -C $OUTPUT/$BASENAME-covRefs.txt $FILEMAME"
        else
            $CMD="$CMD -o $OUTPUT/$BASENAME-transcripts.gtf -A $OUTPUT/$BASENAME-geneAbundances.txt $FILEMAME"
        fi

        printf "$CMD"
        printf "Stringtie assembly for $BASENAME sample finished!!\n"
    done < $1
fi

#printf "Finished all stringtie assemblies!\n\n#####STEP 2####\nRunning merge of all individual assemblies with Cuffmerge..\n"
#cuffmerge -g $REFERENCE_ANNOTATION -s $REFERENCE_GENOME -p 10 -o mergedAssembly individualAssemblies.txt &> log-cuffmerge.txt
#cuffmerge -g $REFERENCE_ANNOTATION  -p 10 -o mergedAssembly individualAssemblies.txt &> log-cuffmerge.txt

printf "Done.\n"

