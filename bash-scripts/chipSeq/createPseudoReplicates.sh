#!/bin/bash

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
	printf "Usage:\n
			1st input tagAlign file	
			2nd output directory
			3rd prefix name for outputFiles.\n"
	exit 1
fi

fileName=$1 # input tagAlign file name
outputDir=$2 # output directory for pseudoreplicate files
outputStub=$3 # prefix name for pseudoReplicate files

nlines=$( zcat ${fileName} | wc -l ) # Number of reads in the tagAlign file
nlines=$(( (nlines + 1) / 2 )) # half that number
zcat "${fileName}" | shuf | split -d -l ${nlines} - "${outputDir}/${outputStub}" # This will shuffle the lines in the file and split it into two parts
gzip "${outputDir}/${outputStub}00"
gzip "${outputDir}/${outputStub}01"
mv "${outputDir}/${outputStub}00.gz" "${outputDir}/${outputStub}.pr1.tagAlign.gz"
mv "${outputDir}/${outputStub}01.gz" "${outputDir}/${outputStub}.pr2.tagAlign.gz"

