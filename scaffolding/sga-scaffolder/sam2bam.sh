#!/bin/bash
#Convert multiple sam files to bam format given a list of files as input

if [ -z $1 ]; then
	printf "Please provide a file with each line representing one file to process.\n"
	exit 1
fi

while read line
do

if [ -f "$line" ]; then
	printf "Processing $line file .."
	samtools view -bS "$line" > ${line/.sam/.bam}	
	printf "Done .."
fi
done < $1
printf "Finished processing all the files"
