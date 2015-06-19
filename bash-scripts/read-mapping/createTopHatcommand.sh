#!/bin/bash
display_usage() { 
printf "Script to generate the topHat2 command to run the mapping of multiple libraries in one run.
Usage:
    -1st argument must be the file with the RNA seq reads to map.
    -2nd argument must be the reference indexed database.
    -3rd argument must be the number of threads to use in the process.
    -4th argument is optional. If enabled, it refers to coverage based search for junctions. Useful for more sensitivity. Available options [true|false].
    -5th argument is optional. If enabled, only the read alignments where both reads in a pair are mapped, are reported. Available options [true|false].\n"
} 


#--no-discordant	For paired reads, report only concordant mappings. Check this options. Only pairs that rely on the insert size provided.
exec="/opt/tools/tophat-2.0.14.Linux_x86_64/tophat2"

#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
	printf "Please provide the arguments required for the script.\n\n"
	display_usage
	exit 1	
fi

index_database="$2"
threads="-p $3"
numb_samples=0
first_pair=true
second_pair=false
pair1=""
pair2=""

#Create command for paired files
while read line
do	

	filename=$line
	if [ -f "$filename" -a "$first_pair" = "true"  ]; then
		pair1="$pair1,$filename"
		first_pair=false
		second_pair=true

	elif [ -f "$filename" -a  "$second_pair" = "true" ]; then
		pair2="$pair2,$filename"
		first_pair=true
		second_pair=false     
		let "numb_samples += 1"	    	

	elif [ ! -f "$filename" ]; then
		echo "ERROR: $filename is not a file.\n"
		display_usage
		exit 1
	fi


done <$1


#Remove first comma in the pairs
pair1=$(echo "$pair1" | sed 's/,/ /')
pair2=$(echo "$pair2" | sed 's/,/ /')

#Finish command
command="$exec $index_database $pair1 $pair2 $threads --b2-very-sensitive --fusion-search"
if [ -n "$4" ]; then
    if [ "$4" = "true" ]  ; then
        command="$command --coverage-search"
    elif [ "$4" != "false" ] ; then
        printf "ERROR: Please provide a valid value for the 4th argument"
    fi
fi

if [ -n "$5" ]; then
    if [ "$5" = "true" ]  ; then
        command="$command --no-mixed"
    elif [ "$5" != "false" ] ; then
        printf "ERROR: Please provide a valid value for the 5th argument"
    fi
fi

#check if output file already exists. If so, delete it
file=run_topHat2_${numb_samples}samples.sh

if [ -f "$file" ]; then
	rm $file
fi
	
##Pass command to script and give permissions to run
echo -e "#!/bin/bash \n$command" > ./$file | chmod +x ./$file






