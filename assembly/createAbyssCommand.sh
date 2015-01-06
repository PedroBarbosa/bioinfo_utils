#!/bin/bash

display_usage() { 
printf "First argument must bee the file. In this file the Paired end libraries need to be first
Second argument must be a flag true/false to use paired end reads to generate the command
Third argument must be a flag true/false to use mate pair reads to generate the command
Fourth argument must be the value of k
fifth argument must be the number of processors to use n
sixth argument is optional. If used, it refers to the minimum number of pairs in different unitigs to be joined into contigs (n)
seventh argument is optional. If used, it refers to the aligner to use (aligner parameter)\n"
} 


exec="abyss-pe"

#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ]; then
	printf "Please provide the arguments required for the script.\n\n"
	display_usage
	exit 1	
else
	k="k=$4"
	np="np=$5"
	name="name=abyss"
	lib="lib='"
	mp="mp='"
	string_samples=""
fi

#check if the optional argument of minimum number of pairs in different unitigs to be joined into contigs is provided (n) [default=10]
if [ -z "$6" ]; then	
	command="$exec $np $k"
else
	min_numb_pairs="n=$6"
	command="$exec $np $k $min_numb_pairs"
fi

#check if the optional argument of the aligner to be used is provided (aligner) [default=map]
if [ -z "$7" ]; then
	echo "No aligner provided. Map will be used as default.."
else
	aligner="aligner=$7"
fi

#add name and aligner (if exists)
command="$command $name $aligner"
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
		let "numb_samples += 1"	  
		lib="$lib lib${numb_samples}"
		string_samples="$string_samples lib${numb_samples}='$pair1 $pair2' "

		  	

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
		let "numb_samples += 1"
		mp="$mp mp${numb_samples}"
		string_samples="$string_samples mp${numb_samples}='$pair1 $pair2'"
	
	
	elif [ ! -f "$filename" -a "$matepairFlag" = "true" ]; then
		echo "$filename" is not a file
		continue
	fi

done <$1

#add last ' to lib and mp and finish the command
if [[ "$2" = "true" && "$3" = "true" ]]; then
	lib="$lib'"
	mp="$mp'"
	command="$command $lib $mp$string_samples"
	
elif [[ "$2" = "true" && "$3" = "false" ]]; then
	lib="$lib'"
	echo $lib
	command="$command $lib$string_samples"
	
elif [[ "$2" = "false" ]]; then	
	echo "Are you trying to assemble without paired end reads ? Please set true for the use of paired end information"

else
	echo "Please provide valid values for the type of libraries to be used"
fi


#echo $command

#check if output file already exists. If so, delete it
file=command_${numb_samples}_samples.sh
echo $file
if [ -f "$file" ]; then
	rm $file
fi
	
echo $file

##Pass command to script and give permissions to run
echo -e "#!/bin/bash \n$command" > ./$file | chmod 755 ./$file

