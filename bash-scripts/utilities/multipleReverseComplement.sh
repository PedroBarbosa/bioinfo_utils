#!/bin/bash

if [ -z "$1" ]; then
	printf "Please provide a file that lists the FASTQ files to process\nSEQTK tool must be on the system path.\nOutput files will be written in the current directory.\n"
	exit 1
fi

pair1=true
pair2=false
#seqkt="/opt/tools/seqtk/seqtk"
while read line
do
	basename=$(basename $line)
	printf "Processing $basename file ..\n"
	if [ "$pair1" == "true" ]; then
		pair1=false
		pair2=true
		if [[ "$basename" == *".fq"* ]] ; then
			output_file=${basename/_1.fq/_RC_1.fq}
		elif [[ "$basename" == *".fastq"* ]] ; then
			output_file=${basename/_1.fastq/_RC_1.fastq}
		else
			printf "Are the files in fastq format? $basename file does not end with .fq or .fastq extension!\n"
			exit 1
		fi
		exec="seqtk seq -r $line"
		$exec > $PWD/$output_file
		
	elif [ "$pair2" == "true" ] ; then
		pair1=true
		pair2=false
                if [[ "$basename" == *".fq"* ]] ; then
                        output_file=${basename/_2.fq/_RC_2.fq}
                elif [[ "$basename" == *".fastq"* ]] ; then
                        output_file=${basename/_2.fastq/_RC_2.fastq}
                else
                        printf "Are the files in fastq format? $basename file does not end with .fq or .fastq extension!\n"
                        exit 1
                fi
		exec="seqtk seq -r $line"
		
		$exec > $PWD/$output_file
	fi
done < $1
