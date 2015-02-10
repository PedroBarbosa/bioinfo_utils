#!/bin/bash

display_usage(){
 printf "Script to run the SGA-align script for multiple libraries, that under the hood runs bwa aln | sampe to generate the alignments'.\n
 Usage:
    -1st argument must be the file. In this file the Paired end libraries need to be first.
    -2nd argument must be a flag to use paired end reads to generate the command.Available option: [true|false].
    -3rd argument must be a flag to use mate pair reads to generate the command.Available option: [true|false].
    -4th argument must be the contigs/scaffolds file in fasta format.
    -5th argument must be the number of threads to use [INT value].
    -6th argument must be the prefix of the output.\n"
}

#################Check if required arguments were provided########
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ]  ; then
    printf "Please provide the required arguments for the script.\n\n"
    display_usage
    exit 1
fi

####################VARIABLES####################
SGA_ALIGN_PATH="/opt/tools/sga/src/bin/sga-align"
numb_samples=0
first_pair=true
second_pair=false
matepairFlag=false


while read line
do

	filename=$line
	#check if reached the MATE PAIR samples, and if not supposed to assemble them, break the loop
	if [[ "$filename" == *"MP"* && "$3" = "true" ]];then
		matepairFlag=true

	elif [[ "$filename" == *"MP"* && "$3" = "false" ]];then
		break
	fi


	if [ -f "$filename" -a "$first_pair" = "true" -a "$2" = "true" -a "$matepairFlag" = "false" ]; then
		pair1=$filename
		first_pair=false
		second_pair=true

	elif [ -f "$filename" -a  "$second_pair" = "true" -a "$2" = "true" -a "$matepairFlag" = "false" ]; then
		pair2=$filename
		first_pair="true"
		second_pair="false"
        printf "Running SGA-align for ${numb_samples} paired end sample ..\n"
		command="$SGA_ALIGN_PATH --name $6${numb_samples} -t $5 $4 $pair1 $pair2"
		$command 2> error.log
		let "numb_samples += 1"

	elif [ ! -f "$filename" -a "$matepairFlag" = "false"  ]; then
		echo "$filename" is not a file
		continue
	fi


	#add mate pair libraries to the command
	if [ "$matepairFlag" = "true" -a -f "$filename" -a "$first_pair" = "true" -a "$3" = "true" ]; then
		pair1=$filename
		first_pair=false
		second_pair=true

	elif [ "$matepairFlag" = "true" -a -f "$filename" -a  "$second_pair" = "true" -a "$3" = "true" ]; then
		pair2=$filename
		first_pair="true"
		second_pair="false"
        printf "Running SGA-align for ${numb_samples} mate pair sample ..\n"
		command="$SGA_ALIGN_PATH --name $6${numb_samples} -t $5 $4 $pair1 $pair2"
		$command 2> error.log
		let "numb_samples += 1"


	elif [ ! -f "$filename" -a "$matepairFlag" = "true" ]; then
		echo "$filename" is not a file
		continue
	fi

done <$1
