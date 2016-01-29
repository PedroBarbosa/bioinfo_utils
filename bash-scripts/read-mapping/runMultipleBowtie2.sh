#!/bin/bash
#Script to run bowtie multiple times for each library. One sam per library will be produced. 
display_usage() { 
printf "First argument must be the file. In this file the Paired end libraries need to be first
Second argument must be a flag true/false to use paired end reads to generate the command
Third argument must be a flag true/false to use mate pair reads to generate the command
Fourth argument must be the reference indexed database
Fifth argument must be the number of threads to use
Sixth argument must be the base name file to write SAM alignments for different libraries
Seventh argument is optional. It refers to the number mismatches allowed in the seed alignments. [Default:0]\n"
} 


exec="bowtie2"

#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ]; then
	printf "Please provide the arguments required for the script.\n\n"
	display_usage
	exit 1	
fi

index_database="-x $4"
threads="-p $5"
sam_basename="-S $6"
numb_samples=0
first_pair=true
second_pair=false
matepairFlag=false


#check if the optional argument of seed coverage is provided
if [ -z "$7" ]; then	
	base_command="$exec $index_database $threads"
else
	mismatches=" -N $7"
	base_command="$exec $index_database $threads $mismatches"
fi

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

	
	#echo $filename
	#add paired end pairs to the command
	if [ -f "$filename" -a "$first_pair" = "true" -a "$2" = "true" -a "$matepairFlag" = "false" ]; then
		pair1="$filename"
		first_pair=false
		second_pair=true

	elif [ -f "$filename" -a  "$second_pair" = "true" -a "$2" = "true" -a "$matepairFlag" = "false" ]; then
		pair2="$filename"
		first_pair=true
		second_pair=false     
		let "numb_samples += 1"
		printf "Running bowtie2 aligner for library ${numb_samples}.."
		sam_file="${sam_basename}_${numb_samples}.sam"
		command="$base_command -1 $pair1 -2 $pair2 $sam_file"	    	
		$command 2> "./stderr.txt"

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
		printf "Running bowtie2 aligner for library ${numb_samples}.."
		sam_file="${sam_basename}_${numb_samples}.sam"
		command="$base_command -1 $pair1 -2 $pair2 $sam_file"
		$command 2> "./stderr.txt"

	
	elif [ ! -f "$filename" -a "$matepairFlag" = "true" ]; then
		echo "$filename" is not a file
		continue
	fi

done <$1
