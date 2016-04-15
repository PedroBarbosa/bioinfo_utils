#!/bin/bash
#Script to generate the bowtie2 command to run the mapping of multiple libraries in one run.
display_usage() { 
printf "First argument must be the file. In this file the Paired end libraries need to be first
Second argument must be a flag true/false to use paired end reads to generate the command
Third argument must be a flag true/false to use mate pair reads to generate the command
Fourth argument must be the reference indexed database
Fifth argument must be the number of threads to use
Sixth argument must be the file to write SAM alignments
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
sam_file="-S $6"
numb_samples=0
first_pair=true
second_pair=false
matepairFlag=false
pair1="-1"
pair2="-2"

#check if the optional argument of seed coverage is provided
if [ -z "$7" ]; then	
	command="$exec $index_database $sam_file $threads"
else
	mismatches=" -N $7"
	command="$exec $index_database $sam_file $threads $mismatches"
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
		pair1="$pair1,$filename"
		first_pair=false
		second_pair=true

	elif [ -f "$filename" -a  "$second_pair" = "true" -a "$2" = "true" -a "$matepairFlag" = "false" ]; then
		pair2="$pair2,$filename"
		first_pair=true
		second_pair=false     
		let "numb_samples += 1"	    	

	elif [ ! -f "$filename" -a "$matepairFlag" = "false"  ]; then
		echo "$filename" is not a file
		continue
	fi


	#add mate pair libraries to the command
	if [ "$matepairFlag" = "true" -a -f "$filename" -a "$first_pair" = "true" -a "$3" = "true" ]; then
		pair1="$pair1,$filename"
		first_pair=false
		second_pair=true

	elif [ "$matepairFlag" = "true" -a -f "$filename" -a  "$second_pair" = "true" -a "$3" = "true" ]; then
		pair2="$pair2,$filename"
		first_pair=true
		second_pair=false      
		let "numb_samples += 1"	    	
	
	
	elif [ ! -f "$filename" -a "$matepairFlag" = "true" ]; then
		echo "$filename" is not a file
		continue
	fi

done <$1


#Remove first comma in the pairs
pair1=$(echo "$pair1" | sed 's/,/ /')
pair2=$(echo "$pair2" | sed 's/,/ /')

#Finish command
command="$command $pair1 $pair2"
#echo $command

#check if output file already exists. If so, delete it
file=command_${numb_samples}_samples.sh

if [ -f "$file" ]; then
	rm $file
fi
	
##Pass command to script and give permissions to run
echo -e "#!/bin/bash \n$command" > ./$file | chmod +x ./$file





