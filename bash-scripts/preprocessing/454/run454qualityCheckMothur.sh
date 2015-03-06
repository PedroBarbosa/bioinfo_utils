#!/bin/bash

display_usage(){
	printf "Script to perform the quality trimming of reads from the Roche 454 platform using Mothur trim.seqs program under the hood. It is required that Mothur is in the system path.

    -1st argument must be a file with the list of FASTA files to process.
    -2nd argument must be a file with the qualities associated. Each line of the 1st and 2nd aragument must represent the same sample.
    -3rd argument must be the number of processors to use [INT].
    -4th argument must be the minimum average quality per read after trimming [0< INT <40].
    -5th argument must be the minimum average quality over a window of size [5th argument] moved in a 1 base step range [INT].
    -6th argument must be the size of the window [INT].
    -7rd argument must be the minimum length of reads to keep after trimming [INT | - ], where '-' means this option will be disabled in Mothur.
    -8th argument must be the maximum length of reads to keep after trimming [INT | - ], where '-' means this option will be disabled in Mothur.
    -9th argument must be the quality threshold for a base call [0< INT <40 | - ] where '-' means this option will be disabled in Mothur. Base calls below this threshold the read will be automatically discarded.
    -10th argument must be the maximum number of ambiguous calls (Ns) allowed in a read [INT | - ] where '-' means this option will be disabled in Mothur.
    -11th argument must be the maximum homopolymer length allowed [INT | - ] where '-' means this option will be disabled in Mothur.\n\n"
}


if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ] || [ -z "$7" ] || [ -z "$8" ] || [ -z "$9" ] || [ -z "${10}" ] || [ -z "${11}" ];then
	printf "Too few arguments.\n"
	display_usage
	exit 1
fi


FASTA_FILE="fasta= "
QUAL_FILE="qfile= "
PROCESSORS="processors= $3"
AVG_READ_QUAL="qaverage= $4"
AVG_WINDOW_QUAL="qwindowaverage= $5"
WINDOW_SIZE="qwindowsize= $6"
EXEC="#trim.seqs($FASTA_FILE, $QUAL_FILE, $PROCESSORS, $AVG_READ_QUAL, $AVG_WINDOW_QUAL, $WINDOW_SIZE"
if [ "$7"  != "-" ]; then
    MIN_LENGTH="minlength= $7"
    EXEC="$EXEC, $MIN_LENGTH"
fi
if [ "$8" != "-" ]; then
    MAX_LENGTH="maxlength= $8"
    EXEC="$EXEC, $MAX_LENGTH"
fi
if [ "$9" != "-" ]; then
    QTHRESHOLD= "qthreshold= $9"
    EXEC="$EXEC, $QTHRESHOLD"
fi
if [ "${10}" != "-" ]; then
    AMBIG_CALLS= "maxambig= ${10}"
    EXEC="$EXEC, $AMBIG_CALLS"
fi
if [ "${11}" != "-" ]; then
    MAX_HOMOPOLYMER= "maxhomop= ${11}"
    EXEC="$EXEC, $MAX_HOMOPOLYMER"
fi

#Close command
EXEC="$EXEC)"

#TRANSFORM TWO FILES IN LISTS AND THEN ACCESS TO THE SAME POSITION OF LISTS TO INSERT AND RUN  THE MOTHUR COMMAND
#ALSO PERFORM SUMMARY.SEQS.

#sort list of sorted and indexed files (increasing insert size)
LIST_FASTA=( $(find $2 -name '*sorted*.bam' | sort -n -t _ -k 2) )
for i in "${list_sorted[@]}"
do
	ORIENTATION="$ORIENTATION $3"
done

