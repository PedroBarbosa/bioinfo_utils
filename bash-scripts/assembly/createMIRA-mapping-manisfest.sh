#!/bin/bash
display_usage(){
 printf "Script to automatically generate the manifest file to input in Mira for mapping assembly projectts.\n
 Usage:
    -1st argument must be the project name to use in Mira.The name will be the prefix for the manifest file generated.
    -2nd argument must be a file representing the reference sequence to map onto.
    -3rd argument must be the name of the reference strain.
    -4rd argument must be a file in which each line represents the path where the sequences of each strain are located with the name of the strain included \
in a tab delimited format.
    -5th argument must be the type of data that will be mapped. Available options: [solexa|454|iontor|sanger]. If solexa (Illumina) paired-end reads are used,\
the file names should include the "_1" and "_2" characters to let the script know the pair info.
    -6rd argument must be the number of threads to use [INT value].
    -7th argument must be the percentage of memory to keep free [0< INT <100}.\n\n

Example of 4th argument:
/directory/where/files/are/located1 strain_x
/directory/where/files/are/located2 strain_y
...
"
}


#################Check if required arguments were provided########
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ] || [ -z "$7" ] ; then
    printf "Please provide the required arguments for the script.\n\n"
    display_usage
    exit 1
fi
####################VARIABLES####################
#Project name
PROJECT_NAME="$1"

#Output file
OUTPUT_FILE="${PROJECT_NAME}-manifest.txt"

#Reference sequence
if [ -f "$2" ]; then
    REF_GENOME="$2"
    REF_STRAIN="$3"
else
    printf "Please providee a valid file representing the reference sequence.\n"
    display_usage
    exit 1

#Sequences to map
TYPE_OF_DATA="$5"
while read line
do
    DATA_PATHS=()
    STRAINS=()
    dir=$(cut -f1 $line)
    strain=$(cut -f2 $line)

    if [ -d "$dir" ] ; then
        DATA_PATHS+=("$dir")
        STRAINS+=("$strain")

    else
        printf "The path $dir does not exist. Please check this out.\n"
        display_usage
        exit 1
    fi
done < "$4"

#Threads
THREADS="$6"
#Memory arguments
MEM_FREE_PERCENT="$7"

###################DEFINING READ GROUPS#############################

function read_groups(){

cat <<EOF >> $OUTPUT_FILE

#Reference genome
readgroup
is_reference
data = $REF_GENOME
strain = $REF_STRAIN
EOF

}

###Data to map########
if [ "$TYPE_OF_DATA" = "solexa" ]; then

    for i in "${!DATA_PATHS[@]}";
    do
    pair1=$(ls ${DATA_PATHS[$i]} | grep "_1")
    pair2=$(ls ${DATA_PATHS[$i]} | grep "_2")
    if [ -z "$pair1" ] || [ -z "$pair2" ]; then
        printf "Please check if the illumina files have the substrings '_1' and/or '_2' in the file names.\n\n"
        display_usage
        exit 1
    else
    cat <<EOF >> $OUTPUT_FILE

#Read group "$i"
readgroup = data_${STRAINS[$i]}
data = ${DATA_PATHS[$i]}$pair1 ${DATA_PATHS[$i]}$pair2
autopairing
segment_placement = ---> <---
technology = $TYPE_OF_DATA
EOF
    fi
    done

else

    for i in "${!DATA_PATHS[@]}";
    do

    cat <<EOF >> $OUTPUT_FILE

#Read group "$i"
readgroup = data_${STRAINS[$i]}
data = ${DATA_PATHS[$i]}
technology = $TYPE_OF_DATA
EOF

    done
fi




function settings(){


cat <<EOF >> $OUTPUT_FILE

#PARAMETERS
parameters = COMMON_SETTINGS -GE:not=$1:amm=no:kpmf=$2 -NW:cmrnl=warn -SK:not=$1:mmhr=4 \\
EOF



}




##################################################################
###################Create command#################################
##################################################################

##START APPENDING GENERAL SETTINGS TO OUTPUT FILE
if [ -f "$OUTPUT_FILE" ] ; then
		rm $OUTPUT_FILE
fi

cat <<EOF >> $OUTPUT_FILE
#Name to the mapping assembly and purpose of the job
project = ${PROJECT_NAME}
job = genome,mapping,accurate
EOF


####READ GROUPS
if [ "$5" = "solexa" ] || [ "$5" = "454" ] || [ "$5" = "iontor" ] || [ "$5" = "sanger" ]; then
    read_groups
else
    printf "Please provide a valid value for the type of data being used in the mapping process.\n"
    display_usage
    exit 1
fi

###PARAMETERS
settings $6 $7
