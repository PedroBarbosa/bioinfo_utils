#!/bin/bash
#This script creates multiple .de files from bam.

display_usage(){
printf "First argument must be the file containing a list of bam files to process.
Second argument is optional. It refers to the minimum number of pairs to make a contig-contig link. Default value is 5.\n"
}

if [ -z "$1" ]; then
	printf "Please provide the required file.\n"
	display_usage
	exit 1
fi

if [ -z "$2" ]; then
	min_pairs=5
else
	min_pairs="$2"
fi

while read line
do	
	printf "Processing $line file ..\n"
	command="/opt/tools/sga/bin/sga-bam2de.pl --prefix ${line/.bam/} -n $min_pairs -m 200 $line"
	$command
done < $1
printf "Done.\n"
