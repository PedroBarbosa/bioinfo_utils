#!/bin/bash

display_usage() { 
printf "First argument must be the file. In this file the Paired end libraries need to be first.
Second argument must be a flag true/false to use paired end reads to generate the command.
Third argument must be a flag true/false to use mate pair reads to generate the command.
Fourth argument must be the scaffolds file.
Fifth argument must be the number of threads to use.
Sixth argument is optional. It refers to the prefix output. Default:[output]
Seventh argument is optional. It refers to the orientation of the mate pair libraries. Available option: [rf|fr]. Default:[rf]\n\n"
} 



function gapclose(){

#Paired end and/or MP reads [mp reads might be treated as IP if they were reversed complemented]. Case letter on IP because it is two different files for each pair
IP=""
OP=""
numb_samples=0
matepairFlag=false
first_pair=true
second_pair=false

#mate pair orientatin provided?
if [ -z "$4" ] || [ "$4" = "rf" ] ; then 
	orientation="outward"
elif [ "$4" = "fr" ]; then
	orientation="inward"
else
	printf "\nPlease provide a correct value for the orientation of the mate pair libraries.\n\n"
	display_usage
	exit 1
fi

while read line
do	
	filename=$line

	#check if reached the MATE PAIR samples, and if not supposed to assemble them, break the loop
	if [[ "$filename" == *"MP"* && "$3" = "true" ]];then
		matepairFlag=true

	elif [[ "$filename" == *"MP"* && "$3" = "false" ]];then
		break
	fi
	
	#paired end
	if [ -f "$filename" -a "$first_pair" = "true" -a "$2" = "true" -a "$matepairFlag" = "false" ]; then
		pair1=$filename
		first_pair=false
		second_pair=true

	elif [ -f "$filename" -a  "$second_pair" = "true" -a "$2" = "true" -a "$matepairFlag" = "false"  ]; then
		pair2=$filename
		first_pair="true"
		second_pair="false"       
		let "numb_samples += 1"	  
		IP="$IP -IP${numb_samples} $pair1 $pair2"

	elif [ ! -f "$filename" -a "$matepairFlag" = "false"  ]; then
		echo "$filename" is not a file
		continue
	fi
	
	#mate pair
	if [ -f "$filename" -a "$first_pair" = "true" -a "$3" = "true" -a "$matepairFlag" = "true" ]; then
                pair1=$filename
                first_pair=false
                second_pair=true

   	elif [ -f "$filename" -a  "$second_pair" = "true" -a "$3" = "true" -a "$matepairFlag" = "true" -a "$orientation" = "inward"  ]; then
                pair2=$filename
                first_pair="true"
                second_pair="false"
                let "numb_samples += 1"
                IP="$IP -IP${numb_samples} $pair1 $pair2"
	
	elif [ -f "$filename" -a  "$second_pair" = "true" -a "$3" = "true" -a "$matepairFlag" = "true" -a "$orientation" = "outward"  ]; then
                pair2=$filename
                first_pair="true"
                second_pair="false"
                let "numb_samples += 1"
                OP="$OP -OP${numb_samples} $pair1 $pair2"

        elif [ ! -f "$filename" -a "$matepairFlag" = "false"  ]; then
                echo "$filename" is not a file
                continue

    	fi

done < "$1"

#add samples to the command
if [ "$3" = "true" -a "$orientation" = "inward" ]; then
	command_gapclose="$command_gapclose $IP"

elif [ "$3" = "true" -a "$orientation" = "outward" ]; then
	
	if [ "$2" = "true" ]; then
		command_gapclose="$command_gapclose $IP $OP"
		
	else 
		command_gapclose="$command_gapclose $OP"
	fi
	
elif [ "$3" = "false" ]; then	
	printf "\nAre you trying to scaffolding without mate pair reads ? Please set true for the use of this reads files.\n\n"
	display_usage
	exit 1

else
	printf "\nPlease provide valid values for the type of libraries to be used.\n\n"
	display_usage
	exit 1
fi

}

generate_gapClose_file(){
	file=command_${numb_samples}_samples_gapClose.sh
        if [ -f "$file" ]; then
                rm $file
        fi

        ##Pass command to script and give permissions to run
        echo -e "#!/bin/bash \n$command_gapclose" > ./$file | chmod +x ./$file
}

exec="/mnt/msa/assembly/platanus-assembly/platanus"
#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] ; then
        printf "\nPlease provide the arguments required for the script.\n\n"
        display_usage
        exit 1
else
      
	scaffolds_file="-c $4"
	threads="-t $5"	
	if [ -z "$6" ] ; then
        	output="-o output"
	else
		output="-o "$6""
	fi       
fi


#generate general command
command_gapclose="$exec gap_close $output $scaffolds_file $threads"
matepairFlag=false

##################################################################
###################Create command#################################
##################################################################
gapclose $1 $2 $3 $7
generate_gapClose_file
