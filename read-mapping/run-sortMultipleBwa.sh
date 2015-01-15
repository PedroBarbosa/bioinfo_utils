#!/bin/bash
#Script to run create bwa index, run bwa and create sorted indexed bam files for each library. One sorted bam per library will be produced. Uses the script created by the besst team and needs bwa to be in the path
display_usage() {
printf "First argument must be the file. In this file the Paired end libraries need to be first
Second argument must be a flag true/false to use paired end reads to generate the command
Third argument must be a flag true/false to use mate pair reads to generate the command
Fourth arguent must be the contigs file to map
Fifth argument must be the number of threads to use in bwa
Sixth argument must be the base name file to write SAM alignments for different libraries\n"
}


#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ]; then
        printf "Please provide the arguments required for the script.\n\n"
        display_usage
        exit 1
fi

exec="python /mnt/msa/scaffolding/besst-scaffolding/reads_to_ctg_map.py"
contigs="$4"
threads="--threads $5"
output_basename="$6"
numb_samples=0
first_pair=true
second_pair=false
matepairFlag=false


#Create command for paired files
while read line
do

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
		printf "\nRunning bwa aligner for library ${numb_samples}..\n"
		command="$exec $pair1 $pair2 $contigs ${output_basename}_${numb_samples} $threads"
		$command

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
                printf "\nRunning bwa aligner for library ${numb_samples}..\n"
                command="$exec $pair1 $pair2 $contigs ${output_basename}_${numb_samples} $threads"
                $command
	
	elif [ ! -f "$filename" -a "$matepairFlag" = "true" ]; then
		echo "$filename" is not a file
		continue
	fi

done <$1
