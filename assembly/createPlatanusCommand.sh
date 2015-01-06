#!/bin/bash

display_usage() { 
printf "First argument must bee the file. In this file the Paired end libraries need to be first
Second argument must be a flag true/false to use paired end reads to generate the command
Third argument must be a flag true/false to use mate pair reads to generate the command
Fourth argument must be the number of threads to use
Fifth argument must be the maximum amount of memory allowed.
Sixth argument is optional. If set to yes, only performs the assemble [default: assembly and scaffolding together]
Seventh argument is optional. If set to yes, only performs scaffolding. The prefix needs to be same as the contigs file, normally 'output' [default: assembly and scaffolding together]\n"
} 

function assembly_contigs() {
	numb_samples=0
	while read line
	do	
		filename=$line

		#check if reached the MATE PAIR samples, and if not supposed to assemble them, break the loop
		if [[ "$filename" == *"MP"* && "$3" = "true" ]];then
			matepairFlag=true

		elif [[ "$filename" == *"MP"* && "$3" = "false" ]];then
			break
		fi
		
		#add files to the command [even though they are paired end there is no difference for the contig step
		if [ -f "$filename" -a "$matepairFlag" = "false" ]; then
			file=$filename
			let "numb_samples += 1"
			string_samples="$string_samples $file"

		elif [ -f "$filename" -a "$matepairFlag" = "true" ]; then
			file=$filename
			let "numb_samples += 1"
			string_samples="$string_samples $file"
		fi
		
	
	done <"$1"
	numb_samples=$(($numb_samples / 2))
	#add samples to the command
	if [[ "$2" = "true" ]]; then
		command_contigs="$command_contigs $string_samples"
	elif [[ "$2" = "false" ]]; then	
		echo "Are you trying to assemble without paired end reads ? Please set true for the use of this read files"
		exit 1
	else
		echo "Please provide valid values for the type of libraries to be used"
		exit 1
	fi

}

function scaffolding(){
contigs_file="-c output_contig.fa"
bubbles_file="-b output_contigBubble.fa"
#Paired end and/or MP reads [mp reads are treated as IP because they are already with the correct orientation]. Case letter on IP because it is two different files for each pair
IP=" "
numb_samples=0
first_pair=true
second_pair=false


while read line
do	
	filename=$line

	#check if reached the MATE PAIR samples, and if not supposed to assemble them, break the loop
	if [[ "$filename" == *"MP"* && "$3" = "true" ]];then
		matepairFlag=true

	elif [[ "$filename" == *"MP"* && "$3" = "false" ]];then
		break
	fi

	if [ -f "$filename" -a "$first_pair" = "true" -a "$2" = "true" ]; then
		pair1=$filename
		first_pair=false
		second_pair=true

	elif [ -f "$filename" -a  "$second_pair" = "true" -a "$2" = "true" ]; then
		pair2=$filename
		first_pair="true"
		second_pair="false"       
		let "numb_samples += 1"	  
		IP="$IP -IP${numb_samples} $pair1 $pair2 "
		  	

	elif [ ! -f "$filename" -a "$matepairFlag" = "false"  ]; then
		echo "$filename" is not a file
		continue
	fi
	
done <"$1"

#add samples to the command
if [[ "$2" = "true" ]]; then
	command_scaffolds="$command_scaffolds $contigs_file $bubbles_file $IP"
elif [[ "$2" = "false" ]]; then	
	echo "Are you trying to scaffolding without paired end reads ? Please set true for the use of this reads files"
	exit 1

else
	echo "Please provide valid values for the type of libraries to be used"
	exit 1
fi
}

function generate_final_file_contigs(){
	#check if output file already exists. If so, delete it
	file=command_${numb_samples}_samples.sh
	if [ -f "$file" ]; then
		rm $file
	fi

	##Pass command to script and give permissions to run
	echo -e "#!/bin/bash \n$command_contigs" > ./$file | chmod 755 ./$file
}

function generate_final_file_scaffolds(){
	#check if output file already exists and if it is supposed to append. If so, delete it
	file=command_${numb_samples}_samples.sh
	if [ -f "$file" ] && [ -z $1 ] ; then #$1 refers to the first parameter passed to the function [in this case, if only scaffolding]. This case in supposed to append, because no information was passed in the script arguments
		echo -e "&& $command_scaffolds" >> ./$file
	elif [[ -f "$file" &&  $1 = "true" ]] ; then
		rm $file
		echo -e "#!/bin/bash \n$command_scaffolds" > ./$file | chmod 755 ./$file
	else
		echo -e "#!/bin/bash \n$command_scaffolds" > ./$file | chmod 755 ./$file
	fi	
}





exec="/mnt/msa/assembly/platanus-assembly/platanus"
#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ]; then
	printf "Please provide the arguments required for the script.\n\n"
	display_usage
	exit 1	
else
	threads="-t $4"
	memory="-m $5"
	output="-o output"
	string_samples="-f"
fi


#generate general command
command_contigs="$exec assemble $output $threads $memory"
command_scaffolds="$exec scaffold $output $threads"
matepairFlag=false



##################################################################
###################Create command#################################
##################################################################
if [ -z "$6" ] && [ -z "$7" ]; then
	assembly_contigs $1 $2 $3
	generate_final_file_contigs 
	scaffolding $1 $2 $3
	generate_final_file_scaffolds $7
elif [ $6 = "yes" ] && [ -z $7 ] ; then
	assembly_contigs $1 $2 $3
	generate_final_file_contigs
elif [ $6 = "yes" ] && [ $7 = "yes" ] ; then
	echo "Please verify whether you want to perform either only assemble or only scaffolding. If you want to perform only scaffolding please set the option for only assemble to 'no' If you want to perform only assemble don't set the 7th argument at all."
	exit 1
elif [ $7 = "yes" ]; then
	scaffolding $1 $2 $3
	generate_final_file_scaffolds $7
else 
	echo "Please check if you provided the right parameters [script is not case sensitive]"
	exit 1
fi

