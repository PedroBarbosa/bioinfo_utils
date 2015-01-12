#!/bin/bash

if [ -z "$1" ]; then
	printf "Please provide the input file.Pairs must come consecutives.\n"
	exit 1
fi
command="sga preprocess --pe-mode 1"
output_file="-o mygenome.fastq"

while read line
do
	command="$command $line"
done < "$1"

#Add output file
command="$command $output_file"

printf "This command is about to run:\n\n$command\n"
$command
