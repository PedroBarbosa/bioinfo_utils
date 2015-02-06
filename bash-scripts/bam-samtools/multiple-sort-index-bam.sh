#!/bin/bash

if [ -z "$1" ]; then
	printf "Please provide the input file containing the bam files to be processsed.\n"
	exit 1
fi	

#memory refers to 100Gb
memory=107374182400
while read line
do	
	printf "Processing $line file..\n"
	output=${line/.sam/-sorted}
	command_view="samtools view -Shub $line" 	
	echo "$command_view" 
	command_sort="samtools sort -m $memory - ${output}"
	echo "$command_sort"
	command_index="samtools index ${output}.bam"
	$command_view | $command_sort 
	printf "Indexing now..\n"
	$command_index
	printf "$line file processed.\n" 
done < "$1"
