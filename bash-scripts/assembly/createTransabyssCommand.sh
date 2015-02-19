#!/bin/bash

display_usage(){
 printf "Script to automatically generate the executable command for the RNA assembly with transabyss. If you want to customize more parameters you can add them manually after this script is ran.\n
 Usage:
    -1st argument must be the project name to use in Transabyss. This will include the name of the output folder and files.
    -2nd argument must be a flag true/false to use paired end reads to generate the command.
    -3rd argument must be a flag true/false to use single end reads to generate the command.
    -4th argument must be the size of the k for the assembly.[INT value]
    -5rd argument must be the number of threads to use [INT value].
    -6th argument must be the file with the paths for the paired end reads to generate the command, if the 2nd argument == true. Pairs must come consecutively. [FILE|'-'], where '-' refers to this argument when no paired reads will be used.
    -7th argument must be the file with the paths for the single end reads to generate the command, if the 3rd argument == true. [FILE|'-'], where '-' refers to this argument when no single end reads will be used.
    -8th argument is optional. It refers to the use of Blat alignments to remove redundant sequences. If set to 'no', no redundancy removal will be performed. [yes|no]. Default:yes.\n\n"
}

#################Check if required arguments were provided########
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ] || [ -z "$7" ] ; then
    printf "Please provide the required arguments for the script.\n\n"
    display_usage
    exit 1
fi

####################VARIABLES####################
#Transabyss exec
TRANSABYSS="/opt/tools/transabyss-master/transabyss"

#Project name
PROJECT_NAME="$1"

#K-mer size
K_SIZE="-k $4"

#Threads to use
THREADS="--threads $5"

#Blat?
if [ -z "$8" ] || [ "$8" == "yes" ]; then
    BLAT='yes'
elif [ "$8" == no ]; then
    BLAT='no'
else
    printf "Please provide a valid valid for the 8th argument, or don't provide it at all, since it is optional.\n\n"
    display_usage
    exit 1
fi

#Single end command
SINGLE_END_DATA="--se"
#Paired end command
PAIRED_END_DATA="--pe"

#Flags
FIRST_PAIR=true
SECOND_PAIR=false



#Single end files to assemble.
if [ "$3" = "true" ] ; then

    if [ -f "$7" ] ; then

        while read line
        do
            if [ -f "$line" ] ; then
                SINGLE_END_DATA="$SINGLE_END_DATA $line"

            else
                printf "File $line does not exist. Please check this out.\n"
                display_usage
                exit 1
            fi
        done < "$7"

    elif [ "$7" = "-" ]; then

        printf "You set the 3rd argument to true, so please provide a file with the single end reads to use.\n"
        display_usage
        exit 1

    else
        printf "Please use a valid value for the 7th argument.\n"
        display_usage
        exit 1
    fi

elif [ "$3" = "false" ] && [ -f "$7" ]; then
    printf "No need to provide a file representing single end reads. 3rd argument was set to false. Please check this out.\n\n"
    display_usage
    exit 1

else
    printf "Please provide valid values for the 3rd argument: [true|false].\n\n"
    display_usage
    exit 1
fi



#Paired end files
if [ "$2" = "true" ] ; then

    if [ -f "$6" ] ; then

        while read line
        do
            if [ -f "$line" ] ; then
                process_paired_end $FIRST_PAIR $SECOND_PAIR $line
            else
                printf "File $line does not exist. Please check this out.\n"
                display_usage
                exit 1
            fi
        done < "$6"

    elif [ "$6" = "-" ]; then

        printf "You set the 2nd argument to true, so please provide a file with the paired end reads to use.\n"
        display_usage
        exit 1

    else
        printf "Please use a valid value for the 6th argument.\n"
        display_usage
        exit 1
    fi

elif [ "$2" = "false" ] && [ -f "$7" ]; then
    printf "No need to provide a file representing paired end reads. 2rd argument was set to false. Please check this out.\n\n"
    display_usage
    exit 1

else
    printf "Please provide valid values for the 2nd argument: [true|false].\n\n"
    display_usage
    exit 1
fi



function process_paired_end(){

    	if [ "$1" = "true" -a "$2" = "false" ]; then
		    pair1=$3
		    FIRST_PAIR="false"
		    SECOND_PAIR="true"
            SINGLE_END_DATA="$SINGLE_END_DATA $pair1"
            
	    elif [ "$1" = "false" -a "$2" = "true" ]; then
		    pair2=$3
		    FIRST_PAIR="true"
		    SECOND_PAIR="false"
		    SINGLE_END_DATA="$SINGLE_END_DATA $pair2"
            PAIRED_END_DATA="$PAIRED_END_DATA $pair1 $pair2"
	    else 
	        printf "Something went wrong with the paired end file processing.\n\n"
	    fi
}




#########Generate final command###########
if [ "$2" = "true" ]; then
    COMMAND="$TRANSABYSS $SINGLE_END_DATA $PAIRED_END_DATA $THREADS $K_SIZE --name $PROJECT_NAME --outdir $PWD/$PROJECT_NAME"
else
    COMMAND="$TRANSABYSS $SINGLE_END_DATA $THREADS $K_SIZE --name $PROJECT_NAME --outdir $PWD/$PROJECT_NAME"
fi

#check if output file already exists. If so, delete it
EXEC_FILE="$PWD/run_transabyss.sh"
if [ -f "$EXEC_FILE" ]; then
	rm $EXEC_FILE
fi

##Pass command to script and give permissions to run
echo -e "#!/bin/bash \nEXEC_FILE" > ./$EXEC_FILE | chmod +x ./$EXEC_FILE