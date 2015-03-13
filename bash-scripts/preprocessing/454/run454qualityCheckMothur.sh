#!/bin/bash

display_usage(){
	printf "Script to perform the quality trimming of reads from the Roche 454 platform using Mothur trim.seqs program under the hood.
It is required that Mothur is in the system path.
Please provide the fasta file with the extension '.fasta', otherwise the script will not work properly.\n

    -1st argument must be a file with the list of FASTA files to process.
    -2nd argument must be a file with the qualities associated. Each line of the 1st and 2nd argument must represent the same sample.
    -3rd argument must be the number of processors to use [INT].
    -4th argument must be the minimum average quality per read after trimming [0< INT <40].
    -5th argument must be the minimum average quality over a window of size [5th argument] moved in a 1 base step range [INT].
    -6th argument must be the size of the window [INT].
    -7rd argument must be the minimum length of reads to keep after trimming [INT | - ], where '-' means this option will be disabled in Mothur.
    -8th argument must be the maximum length of reads to keep after trimming [INT | - ], where '-' means this option will be disabled in Mothur.
    -9th argument must be the quality threshold for a base call [0< INT <40 | - ] where '-' means this option will be disabled in Mothur. Base calls below this threshold will automatically discard the read.
    -10th argument must be the maximum number of ambiguous calls (Ns) allowed in a read [INT | - ] where '-' means this option will be disabled in Mothur.
    -11th argument must be the maximum homopolymer length allowed [INT | - ] where '-' means this option will be disabled in Mothur.\n\n"
}


if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ] || [ -z "$7" ] || [ -z "$8" ] || [ -z "$9" ] || [ -z ${10} ] || [ -z ${11} ];then
	printf "Too few arguments.\n"
	display_usage
	exit 1
fi


FASTA_FILE="fasta="
QUAL_FILE="qfile="
PROCESSORS="processors=$3"
AVG_READ_QUAL="qaverage=$4"
AVG_WINDOW_QUAL="qwindowaverage=$5"
WINDOW_SIZE="qwindowsize=$6"
EXEC="trim.seqs($FASTA_FILE,$QUAL_FILE,$PROCESSORS,$AVG_READ_QUAL,$AVG_WINDOW_QUAL,$WINDOW_SIZE"
if [ "$7"  != "-" ]; then
    MIN_LENGTH="minlength=$7"
    EXEC="$EXEC,$MIN_LENGTH"
fi
if [ "$8" != "-" ]; then
    MAX_LENGTH="maxlength=$8"
    EXEC="$EXEC,$MAX_LENGTH"
fi
if [ "$9" != "-" ]; then
    QTHRESHOLD="qthreshold=$9"
    EXEC="$EXEC,$QTHRESHOLD"
fi
if [ "${10}" != "-" ]; then
    AMBIG_CALLS="maxambig=${10}"
    EXEC="$EXEC,$AMBIG_CALLS"
fi
if [ "${11}" != "-" ]; then
    MAX_HOMOPOLYMER= "maxhomop=${11}"
    EXEC="$EXEC,$MAX_HOMOPOLYMER"
fi

#Add output dir and close command
OUT_DIR="outputdir=$PWD"
EXEC="$EXEC,$OUT_DIR)"



#Transform lines in files into arrays
readarray LIST_FASTA < "$1"
readarray LIST_QUAL < "$2"

#Print parameters set
printf "\nMothur will run with the following settings for all files:
Number of processors:   $3
Minimum average quality for a read: $4
Minimum average quality over a window:    $5
Size of the window: $6\n"
if [ "$7"  == "-" ]; then
    printf "Minimum length for a read:  Disabled\n"
else
    printf "Minimum length for a read:  $7\n"
fi

if [ "$8"  == "-" ]; then
    printf "Maximum length for a read:  Disabled\n"
else
    printf "Maximum length for a read:  $8\n"
fi

if [ "$9" == "-" ]; then
    printf "Quality threshold for a base call:  Disabled\n"
else
    printf "Quality threshold for a base call:  $9\n"
fi

if [ "${10}"  == "-" ]; then
    printf "Maximum number of ambiguous calls allowed (N's): Disabled\n"
else
    printf "Maximum number of ambiguous calls allowed (N's): ${10}\n"
fi

if [ "${11}"  == "-" ]; then
    printf "Maximum length of homopolymer allowed:  Disabled\n\n"
else
    printf "Maximum length of homopolymer allowed:  ${11}\n\n"
fi
printf "Most of the running time will take place in the conversion of the output files to fastq format. This is because 'make.fastq' program\
 within Mothur does not have an option to run in several threads.\n\n"


#Process and run each sample
for i in "${!LIST_FASTA[@]}"
do
    #files to process
	NEW_FASTA_FILE="fasta=${LIST_FASTA[$i]}"
    NEW_QUAL_FILE="qfile=${LIST_QUAL[$i]}"

    #remove white char from string, otherwise mothur will fail
    NEW_FASTA_FILE=$(echo "$NEW_FASTA_FILE" | sed 's/[[:space:]]//g')
    NEW_QUAL_FILE=$(echo "$NEW_QUAL_FILE" | sed 's/[[:space:]]//g')

    #exec command ready
    EXEC="${EXEC/$FASTA_FILE,$QUAL_FILE/$NEW_FASTA_FILE,$NEW_QUAL_FILE}"
    echo $SHORT_NEW_FASTA_FILE
    #remove full path
    SHORT_NEW_FASTA_FILE=$(basename ${NEW_FASTA_FILE})
    SHORT_NEW_QUAL_FILE=$(basename ${NEW_QUAL_FILE})

    echo $SHORT_NEW_FASTA_FILE

    #create the batch file to run mothur programs at once. trim.seqs [quality trimming], summary.seqs [produce trimming points file] and make.fastq [produce final fastq file]
    echo -e "#File to run mothur in batch mode.\n$EXEC\nsummary.seqs()\nsummary.seqs(fasta=${SHORT_NEW_FASTA_FILE/.fasta/.scrap.fasta})\nmake.fastq(\
fasta=${SHORT_NEW_FASTA_FILE/.fasta/.trim.fasta},qfile=${SHORT_NEW_QUAL_FILE/.qual/.trim.qual})" > "$PWD/batch_file.txt"


    #Run mothur
    printf "Running mothur on ${NEW_FASTA_FILE/fasta=/} sample ..\n"
    RUN_EXEC="mothur $PWD/batch_file.txt"
    LOG_BASENAME=${SHORT_NEW_FASTA_FILE/.fasta/.txt}
    $RUN_EXEC &> "$PWD/log_${LOG_BASENAME}"

    #change the name of final fastq output file
    mv "${SHORT_NEW_FASTA_FILE/.fasta/}.trim.fastq" "${SHORT_NEW_FASTA_FILE/.fasta/}.mothur-output.fastq"

    #remove all log files produced by mothur
    rm *.logfile

    #exec command reset
    EXEC="${EXEC/$NEW_FASTA_FILE,$NEW_QUAL_FILE/$FASTA_FILE,$QUAL_FILE}"
    printf "Done.\n\n"
done


