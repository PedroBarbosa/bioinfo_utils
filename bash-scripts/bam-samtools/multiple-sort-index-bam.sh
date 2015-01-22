#!/bin/bash

if [ -z "$1" ]; then
	printf "Please provide the input file containing the bam files to be processsed.\n"
	exit 1
fi	

#memory refers to 100Gb
memory=107374182400
while read line
do	
	printf "Sorting $line file..\n"
	output=${line/.bam/-sorted}	
	command_sort="samtools sort -m $memory $line $output"
	$command_sort 
	printf "Indexing it now..\n"
	command_index="samtools index ${output}.bam"
	$command_index
	printf "$line file processed.\n" 
done < "$1"
