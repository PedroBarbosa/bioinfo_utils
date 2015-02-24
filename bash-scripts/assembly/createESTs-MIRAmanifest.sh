#!/bin/bash
display_usage(){
 printf "Script to automatically generate the manifest file to input in Mira for EST assembly projectts.\n
 Usage:
    -1st argument must be the project name to use in Mira.The name will be the prefix for the manifest file generated.
    -2nd argument must be a file with the list of files to assemble.
    -3rd argument must be the number of threads to use [INT value].
    -4th argument must be the maximum ammount of memory to use [INT value].
    -5th argument must be the percentage of memory to keep free [0< INT <100}.
    -6th argument is optional. If set to yes, Mira will force some steps to use less memory, with the cost in the runtime. [yes|no] Default:no.
    -7th argument is optional. If set to yes, Mira will execute the clipping algorithm on the reads along with the Poly-A tail removal. [yes|no] Default: yes.\n\n"
}
#You should always use the resume option ('-r') when calling Mira. It will resume the assembly at the point where some special files were written.\n

#################Check if required arguments were provided########
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] ; then
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
    if [ -f "$line" ] ; then
        ESTs_454_DATA="$ESTs_454_DATA $line"

    else
        printf "File $line does not exist. Please check this out.\n"
        display_usage
        exit 1
    fi
done < "$2"




###################DEFINING READ GROUPS#############################

function read_groups(){


###454 ESTs ########

cat <<EOF >> $OUTPUT_FILE

#ESTs reads group
readgroup = 454-ESTs-data
data = ${ESTs_454_DATA}
technology = 454
EOF

}


function settings(){


####PARAMETERS###
#Min reads per contigs [default 2]. Should i put one to get very low abundant transcripts? test it later : −AS:mrpc = 1
#Is poly tails removed by seqclean ? If so, disable this in Mira
#Repeat mask parameters. [-HS:mnr] and [-HS:nrr] respectively [-HS:nrc]. I'll come back to [-SK:bph]"
#mmhr=10 megahubs ratio - default 0. For ESTs projects it should be changed otherwise MIRA will crash in the middle of the process

if [ -z "$4" ] || [ "$4" = "no" ]; then
    printf "Not forcing Mira to use less memory.\n"
cat <<EOF >> $OUTPUT_FILE

#PARAMETERS
parameters = COMMON_SETTINGS -GE:not=$1:amm=no:mps=$2:kpmf=$3 -NW:cmrnl=warn -SK:not=$1:mmhr=10 \\
EOF

#add parameters to force memory reduction
else
    print "Forcing Mira to use less memory.\n"
cat <<EOF >> $OUTPUT_FILE

#PARAMETERS
parameters = COMMON_SETTINGS -GE:not=$1:amm=no:mps=$2:kpmf=$3 -NW:cmrnl=warn -SK:not=$1:mmhr=10:mhpr=500:mhim=10000000 \\
EOF
fi

##clipping algorithm [disable both poly-A tail removal step and the proposed end clipping algorithm]
if [ -z "$5" ] || [ "$5" = "yes" ]; then
    printf "Clipping will be performed as default.\n"
else
    printf "Clipping will not be performed.\n"
cat <<EOF >> $OUTPUT_FILE
454_SETTINGS -CL:pec=no:cpat=no

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
if [ -z "$6" ] && [ -z "$7" ]; then
    settings $3 $4 $5
elif [ "$6" = "yes" ] || [ "$6" = "no" ] && [ -z "$7" ]; then
    settings $3 $4 $5 $6
elif [ "$6" = "yes" ] || [ "$6" = "no" ] && [ "$7" = "yes" ] || [ "$7" = "no" ] ; then
    settings $3 $4 $5 $6 $7
else
    printf "\nPlease provide a rigth value for the 6th or 7th parameter.\n\n"
    display_usage
    exit 1
fi

