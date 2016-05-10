#!/bin/bash

display_usage() { 
printf "First argument must bee the file. In this file the Paired end libraries need to be first
Second argument must be a flag true/false to use paired end reads to generate the command
Third argument must be a flag true/false to use mate pair reads to generate the command
Fourth argument must be the value of k
fifth argument must be the number of processors to use n
sixth argument is optional. If used, it refers to the minimum seed coverage (min cov required)\n"
} 


exec="/opt/tools/Ray-2.3.1/Ray"

#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ]; then
	printf "Please provide the arguments required for the script.\n\n"
	display_usage
	exit 1	
else
	k="-k $4"
	n="-n $5"
fi

#check if the optional argument of seed coverage is provided
if [ -z "$6" ]; then	
	command="mpiexec $n $exec $k"
else
	min_coverage="-use-minimum-seed-coverage $6"
	command="mpiexec $n $exec $k $min_coverage"
fi

command="$command -disable-scaffolder -write-marker-summary -show-memory-usage"
numb_samples=0
first_pair=true
second_pair=false
matepairFlag=false


#Create command for paired files
while read line
do	
	#remove leading and trailing spaces
	#filename=$(echo $line | sed -e 's/^ *//' -e 's/ *$//')
	#path=/mnt/msa/workflow_scripts/
	#filename=$path$line
	filename=$line

	#check if reached the MATE PAIR samples, and if not supposed to assemble them, break the loop
	if [[ "$filename" == *"MP"* && "$3" = "true" ]];then
		matepairFlag=true

	elif [[ "$filename" == *"MP"* && "$3" = "false" ]];then
		break
	fi

	
	#echo $filename
	#add paired end pairs to the command
	if [ -f "$filename" -a "$first_pair" = "true" -a "$2" = "true" -a "$matepairFlag" = "false" ]; then
		pair1=$filename
		first_pair=false
		second_pair=true

	elif [ -f "$filename" -a  "$second_pair" = "true" -a "$2" = "true" -a "$matepairFlag" = "false" ]; then
		pair2=$filename
		first_pair="true"
		second_pair="false"       
		command="$command -p $pair1 $pair2"
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
		command="$command -p $pair1 $pair2"
		let "numb_samples += 1"	    	
	
	
	elif [ ! -f "$filename" -a "$matepairFlag" = "true" ]; then
		echo "$filename" is not a file
		continue
	fi

done <$1


#Add output directory
command="$command -o ./output"
#echo $command

#check if output file already exists. If so, delete it
file=command_${numb_samples}_samples.sh

if [ -f "$file" ]; then
	rm $file
fi
	
##Pass command to script and give permissions to run
echo -e "#!/bin/bash \n$command" > ./$file | chmod 755 ./$file

