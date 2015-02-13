#!/bin/bash
"[-HS:mnr] and [-HS:nrr] respectively [-HS:nrc]. I'll come back to [-SK:bph]"
display_usage(){
 printf "Script to automatically generate the manifest file to input in Mira for an EST assembly projectts.\n
 Usage:
    -1st argument must be the project name to use in Mira.The name will be the prefix for the manifest file generated.
    -2nd argument must be a file with the list of files to assemble.
    -3rd argument must be the number of threads to use [INT value].
    -4th argument must be the maximum ammount of memory to use [INT value].
    -5th argument must be the percentage of memory to keep free [0< INT <100}.
    -6th argument is optional. If set to yes, Mira will force some steps to use less memory, with the cost in the runtime. [yes|no] Default:no.\n"
}
#You should always use the resume option ('-r') when calling Mira. It will resume the assembly at the point where some special files were written.\n

#################Check if required arguments were provided########
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ] || [ -z "$7" ] || [ -z "$8" ] ; then
    printf "Please provide the required arguments for the script.\n\n"
    display_usage
    exit 1
fi
####################VARIABLES####################
#Project name
PROJECT_NAME="$1"

#Output file
OUTPUT_FILE="${PROJECT_NAME}-manifest.txt"

#454 files to assemble.
while read line
do
ESTs_454_DATA="$ESTs_454_DATA $line"
done < "$2"




###################DEFINING READ GROUPS#############################

function read_groups(){


###454 ESTs ########
if [ -d "$ESTs_454_DATA" ] ; then
    cat <<EOF >> $OUTPUT_FILE

#ESTs reads group
readgroup = 454-ESTs-data
data = ${ESTs_454_DATA}
technology = 454
EOF

else
    printf "File regarding the path for the 454 fasta files is not valid. Change the '$ESTs_454_DATA' variable in
        the script to a valid path.\n"
            display_usage
            exit 1
fi

}

function settings(){
####PARAMETERS###



#Min reads per contigs [default 2]. Should i put one to get very low abundant transcripts? test it later : −AS:mrpc = 1
#Repeat mask parameters. [-HS:mnr] and [-HS:nrr] respectively [-HS:nrc]. I'll come back to [-SK:bph]"

if [ -z "$4" ] || [ "$4" = "no" ]; then
    printf "Not forcing Mira to use less memory.\n"
cat <<EOF >> $OUTPUT_FILE

#PARAMETERS
parameters = COMMON_SETTINGS -GE:not=$1:amm=no:mps=$2:kpmf=$3 -NW:cmrnl=warn -SK:not=$1 \

EOF

#add parameters to force memory reduction
else
cat <<EOF >> $OUTPUT_FILE

#PARAMETERS
parameters = COMMON_SETTINGS -GE:not=$1:amm=no:mps=$2:kpmf=$3 -NW:cmrnl=warn -SK:not=$1:mhpr=500:mhim=10000000 \

EOF
fi
}




##################################################################
###################Create command#################################
##################################################################

##START APPENDING GENERAL SETTINGS TO OUTPUT FILE
if [ -f "$OUTPUT_FILE" ] ; then
		rm $OUTPUT_FILE
fi

cat <<EOF >> $OUTPUT_FILE
#Name to the assembly and purpose of the job
project = ${PROJECT_NAME}
job = est,denovo,accurate
EOF


####READ GROUPS
read_groups

###PARAMETERS
if [ -z "$6" ]; then
    settings $4 $5 $6
elif [ "$6" = "yes" ] || [ "$6" = "no" ]; then
    settings $6 $7 $8 $6
else
    printf "\nPlease provide a rigth value for the 9th parameter.\n\n"
    display_usage
    exit 1
fi

