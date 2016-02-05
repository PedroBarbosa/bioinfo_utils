#!/bin/bash
#Script to run bwa for each library with default parameters. One sam per library will be produced. 
display_usage() { 
printf "First argument must be the file. In this file the Paired end libraries need to be first
Second argument must be a flag true/false to use paired end reads to generate the command
Third argument must be a flag true/false to use mate pair reads to generate the command
Fourth argument must be the prefix for the reference indexed database
Fifth argument must be the number of threads to use
Sixth argument must be the base name file to write SAM alignments for different libraries [REMOVED]\n"
} 


exec="bwa mem"

#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ]; then
	printf "Please provide the arguments required for the script.\n\n"
	display_usage
	exit 1	
fi

index_database="$4"
threads="-t $5"
#sam_basename=" $"
numb_samples=0
first_pair=true
second_pair=false
matepairFlag=false


base_command="$exec $threads $index_database "

#Create command for paired files
while read line
do	
	#remove leading and trailing spaces
	#filename=$(echo $line | sed -e 's/^ *//' -e 's/ *$//')
	#path=/mnt/msa/workflow_scripts/
	#filename=$path$line
	filename=$line
	
	#check if reached the MATE PAIR samples, and if not supposed to map them, break the loop
	if [[ "$filename" == *"MP"* && "$3" = "true" ]];then
		matepairFlag=true

	elif [[ "$filename" == *"MP"* && "$3" = "false" ]];then
		break
	fi

	
#	echo $filename
	#add paired end pairs to the command
	if [ -f "$filename" -a "$first_pair" = "true" -a "$2" = "true" -a "$matepairFlag" = "false" ]; then
		pair1="$filename"
		first_pair=false
		second_pair=true
	elif [ -f "$filename" -a  "$second_pair" = "true" -a "$2" = "true" -a "$matepairFlag" = "false" ]; then
		pair2="$filename"

		sam_basename=$(basename $pair2 | cut -d "_" -f1)
		first_pair=true
		second_pair=false     
		let "numb_samples += 1"
		printf "Running bwa mem aligner for library ${numb_samples}..\n" 2> "log.txt"
		sam_file="${sam_basename}.sam"
		command="$base_command $pair1 $pair2"
		$command 1> $sam_file 2>"log.txt"
	elif [ ! -f "$filename" -a "$matepairFlag" = "false"  ]; then
		echo "$filename" is not a file
		continue
	fi


	#add mate pair libraries to the command
	if [ "$matepairFlag" = "true" -a -f "$filename" -a "$first_pair" = "true" -a "$3" = "true" ]; then
		pair1="$filename"
		first_pair=false
		second_pair=true

	elif [ "$matepairFlag" = "true" -a -f "$filename" -a  "$second_pair" = "true" -a "$3" = "true" ]; then
		pair2="$filename"
		first_pair=true
		second_pair=false      
		let "numb_samples += 1"	 
		printf "Running bwa mem aligner for library ${numb_samples}..\n" 
		sam_file="${sam_basename}_${numb_samples}.sam"
		command="$base_command $pair1 $pair2" 
		echo $command
		$command 1> $sam_file 2>"log.txt"

	
	elif [ ! -f "$filename" -a "$matepairFlag" = "true" ]; then
		echo "$filename" is not a file
		continue
	fi
done <$1
