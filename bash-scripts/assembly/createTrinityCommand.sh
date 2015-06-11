#!/bin/bash
display_usage(){
 printf "Script to automatically generate the executable command for the RNA seq assembly with Trinity. If you want to customize more parameters you can add them manually after this script is ran.\n
 Usage:
    -1st argument must be the project name to use in Transabyss. This will include the name of the output folder and files.
    -2st argument must a file with the path for the paired end files. Pairs must come consecutively in file.
    -3nd argument must be the number of CPU threads to use [INT value].
    -4rd argument must be the maximum ammount of memory to use [INT value].
    -5th argument is optional. It refers to the read orientation. Availabe options: [RF|FR].
    -6rd argument is optional. It refers to the reads normalization step. Available options [true|false].
    -7th argument is optional. Perform genome_guided assembly when reference genome is available. INPUT: [bam file].\n\n"
}

#################Check if required arguments were provided########
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] ; then
    printf "Please provide the required arguments for the script.\n\n"
    display_usage
    exit 1
fi

####################VARIABLES####################
#Transabyss exec
TRINITY="/opt/tools/trinity-2.0.6/Trinity"

#Project name
PROJECT_NAME="--output $PWD/$1"

#File of paired end reads
FILE="$2"
#Threads to use
THREADS="--CPU $3"

#Max memory to use
MEMORY="--max_memory $4"

#Flags
FIRST_PAIR="--left "
SECOND_PAIR="--right "


function process_paired_end(){
first_pair=true
second_pair=false
numb_samples=0
while read line
do
    filename=$line
    if [ -f "$filename" -a "$first_pair" = "true" ] ; then
        pair1="$filename"
        first_pair="false"
        second_pair="true"
        FIRST_PAIR="$FIRST_PAIR$pair1,"

	elif [ -f "$filename" -a  "$second_pair" = "true" ]; then
		pair2="$filename"
		first_pair="true"
		second_pair="false"
		SECOND_PAIR="$SECOND_PAIR$pair2,"
		let "numb_samples += 1"

	elif [ ! -f "$filename" ]; then
		echo "ERROR: $filename is not a file.\n\n"
		display_usage
		exit 1
	fi

done <$1
printf "$numb_samples pairs will be used to feed trinity\n."
}



#########Generate final command###########
if [ -z "$5" ] ; then
    process_paired_end "$FILE"
    COMMAND="$TRINITY --seqType fq $MEMORY $THREADS $FIRST_PAIR $SECOND_PAIR $PROJECT_NAME"
elif [ -z "$6" ] ; then
    if [ "$5" != "RF" ] || [ "$5" != "FR" ]; then
        printf "ERROR: Please provide a valid value for the orientation of the reads.\n\n"
        display_usage
        exit 1
    else
        COMMAND="$TRINITY --seqType fq $MEMORY $THREADS $FIRST_PAIR $SECOND_PAIR $PROJECT_NAME --SS_lib_type $5"
    fi
elif [ -z "$7" ]; then
    if [ "$6" != "true" ] || [ "$6" != "false" ]; then
        printf "ERROR: Please provide a valid value for the normalization parameter.\n\n"
        display_usage
        exit 1
    elif [ "$6" = "true" ]; then
        COMMAND="$TRINITY --seqType fq $MEMORY $THREADS $FIRST_PAIR $SECOND_PAIR $PROJECT_NAME --SS_lib_type $5 --normalize_reads"
    fi
elif [ -n "$8" ]; then
    COMMAND="$TRINITY --seqType fq $MEMORY $THREADS $FIRST_PAIR $SECOND_PAIR $PROJECT_NAME --SS_lib_type $5 --normalize_reads --genome_guided_bam $8"

fi

#check if output file already exists. If so, delete it
EXEC_FILE="$PWD/run_trinity.sh"
if [ -f "$EXEC_FILE" ]; then
	rm "$EXEC_FILE"
fi

##Pass command to script and give permissions to run
printf "Done.\n"
echo -e "#!/bin/bash \n$COMMAND" > "$EXEC_FILE" | chmod +x "$EXEC_FILE"