#!/bin/bash
display_usage(){
 printf "Script to automatically generate the executable command for the RNA seq assembly with Trinity. If you want to customize more parameters you can add them manually after this script is ran.\n
 Usage:
    -1st argument must be the project name to use in Trinity. This will include the name of the output folder and files.
    -2st argument must a file with the path for the paired end files. Pairs must come consecutively in file.
    -3nd argument must be the number of CPU threads to use [INT value].
    -4rd argument must be the maximum ammount of memory to use [INT value].
    -5th argument is optional. Is this RNA-seq study stranded specific ? If so, what's the read orientation? Availabe options: [RF|FR|unstranded]. Default:'unstranded'.
    -6rd argument is optional. It refers to the reads normalization step. If set to true, reads beyond 75 of coverage will be discaraded. Available options [true|false].
    -7th argument is optional. Perform genome_guided assembly when reference genome is available. If BAM passed, a maximum genome intron size will be set to 15000, as trinity requires this parameter when running on genome guided mode. INPUT: [bam file].\n\n"
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
PROJECT_NAME="--output $PWD/$1-trinity"

#File of paired end reads
FILE="$2"
#Threads to use
THREADS="--CPU $3"

#Max memory to use
MEMORY="--max_memory $4G"

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

#elif [ ! -f "$filename" ]; then
#		echo "ERROR: $filename is not a file.\n\n"
#		display_usage
#		exit 1
#	fi
fi
done <$1
FIRST_PAIR=${FIRST_PAIR::-1}
SECOND_PAIR=${SECOND_PAIR::-1}

}



#########Generate final command###########
if [ -z "$5" ] ; then
    process_paired_end "$FILE"
    COMMAND="$TRINITY --seqType fq --min_kmer_cov 2 $MEMORY $THREADS $FIRST_PAIR $SECOND_PAIR $PROJECT_NAME"

elif [ "$5" = "RF" ] || [ "$5" = "FR" ]; then
    process_paired_end "$FILE"
    COMMAND="$TRINITY --seqType fq --min_kmer_cov 2 $MEMORY $THREADS $FIRST_PAIR $SECOND_PAIR $PROJECT_NAME --SS_lib_type $5"
elif [ "$5" = "unstranded" ]; then
    process_paired_end "$FILE"
    COMMAND="$TRINITY --seqType fq --min_kmer_cov 2 $MEMORY $THREADS $FIRST_PAIR $SECOND_PAIR $PROJECT_NAME"

else
    printf "ERROR: Please provide a valid value for the strand specificity of the libraries.\n\n"
    display_usage
    exit 1
fi

if [ -n "$6" ] ; then
    if [ "$6" == 'true' ] ; then
        COMMAND="$COMMAND --normalize_reads --normalize_max_read_cov 75"

    elif [ "$6" != "false" ]; then
        printf "ERROR: Please provide a valid value for the normalization parameter.\n\n"
        display_usage
        exit 1
    fi
fi


if [ -n "$7" ]; then
    if [ -f "$7" ] ; then
        COMMAND="$TRINITY --genome_guided_bam $7 --genome_guided_max_intron 15000 $MEMORY $THREADS"
	printf "As this is a genome guided assembly, the reads files will not be used. The BAM file is the only requirement.\n"
    else
        printf "ERROR: Please provide a valid file for the genome_guided parameter.\n\n"
        display_usage
        exit 1
    fi

fi

#check if output file already exists. If so, delete it
EXEC_FILE="$PWD/run_trinity.sh"
if [ -f "$EXEC_FILE" ]; then
	rm "$EXEC_FILE"
fi

##Pass command to script and give permissions to run
printf "Command generated.\n"
printf "If not a genome guided assembly,$numb_samples pairs of reads will be used to feed trinity\n"
echo -e "#!/bin/bash \n$COMMAND" > "$EXEC_FILE" | chmod +x "$EXEC_FILE"

