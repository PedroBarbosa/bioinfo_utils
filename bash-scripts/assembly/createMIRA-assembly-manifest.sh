#!/bin/bash

display_usage(){
 printf "Script to automatically generate the manifest file to input in Mira for de novo assembly projects. Paired end reads\
 must be in the traditional 'forward-reverse' orientation while mate pair must come in the 'reverse-forward' orientation'.\n
 Usage:
    -1st argument must be the project name to use in Mira.The name will be the prefix for the manifest file generated.
    -2nd argument must be a flag to use paired end reads to generate the command.Available option: [true|false].
    -3rd argument must be a flag to use mate pair reads to generate the command.Available option: [true|false].
    -4th argument must be a flag to use 454 reads reads to generate the command.Available option: [true|false].
    -5th argument must be a flag to Mira auto estimate insert sizes of libraries: [true|false].
    -6th argument must be the number of threads to use [INT value].
    -7th argument must be the maximum ammount of memory to use [INT value].
    -8th argument must be the percentage of memory to keep free [0< INT <100].
    -9th argument is optional. If set to true, Mira will force some steps to use less memory, with the cost in the runtime. [true|false] Default:false.
    -10th argument is optional. If set to true, Mira will use IonTorrent reads to perform assembly. [true|false] Default:false.\n"

}
#You should always use the resume option ('-r') when calling Mira. It will resume the assembly at the point where some special files were written.\n

#################Check if required arguments were provided########
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ] || [ -z "$7" ] || [ -z "$8" ] ; then
    printf "ERROR: Please provide the required arguments for the script.\n\n"
    display_usage
    exit 1
fi
####################VARIABLES####################
#Project name
PROJECT_NAME="$1"

#Output file
OUTPUT_FILE="${PROJECT_NAME}-manifest.txt"

#Path to the genomic 454 sequencing FASTA files - All files of the directory will be used so we can provide a path in Mira
GENOMIC_454_PATH="/mnt/msa/BIOCANT/genomic-data/SFF_genom/FASTQ_RAW/"

#Path to the genomic IonTOrrent sequencing FASTQ file
ION_TORRENT_PATH="/mnt/msa/celia_Leao_INIAV/old-analysis/Raw-data/strain-3A"

#Path to the file that lists Illumina files wit the threshold Q20L80/Q20L40
#LIST_ILLUM_PE_MP_PATH="/mnt/msa/workflow_scripts/LIST_FILES/listFiles_Q20L80-PE_Q20L20-MP.txt"
LIST_ILLUM_PE_MP_PATH="/mnt/msa/workflow_scripts/LIST_FILES/listFiles_Q35L90-PE_Q35L40-MP.txt"

##Fragment sizes estimations for each type of library
PE170_INSZ_MIN_MAX="100 450 exclusion_criterion autorefine"
PE500_INSZ_MIN_MAX="350 800 exclusion_criterion autorefine"
PE800_INSZ_MIN_MAX="650 1100 exclusion_criterion autorefine"
MP2000_INSZ_MIN_MAX="1750 2400 exclusion_criterion autorefine"
MP5000_INSZ_MIN_MAX="4700 5300 exclusion_criterion autorefine"

##Data for each read group
PE170=""
PE500=""
PE800=""
MP2000=""
MP5000=""




###################DEFINING READ GROUPS#############################

function read_groups(){

#If assemble paired end
if [ "$2" = "true" ]; then
    #check if is file
 
    if [ -f ${LIST_ILLUM_PE_MP_PATH} ]; then

        #loop to separate PE libraries by insert size to create the read group
        while read line
        
        do
          
            if [[ "$line" == *"PE170"* ]] ; then
                PE170="${PE170}fastq::$line ";
            elif [[ "$line" == *"PE500"* ]]; then
                PE500="${PE500}fastq::$line ";
            elif [[ "$line" == *"PE800"* ]] ; then
                PE800="${PE800}fastq::$line ";

            else
                continue
            fi
        done < ${LIST_ILLUM_PE_MP_PATH}

        #check the auto estimate parameter
        if [ "$5" = "true" ] ; then
            cat <<EOF >> $OUTPUT_FILE
            
#pe170 read group
readgroup = illumina-pe-170
data = ${PE170}
autopairing
segment_placement = ---> <---
technology = solexa

#pe500 read group
readgroup = illumina-pe-500
data = ${PE500}
autopairing
segment_placement = ---> <---
technology = solexa

#pe800 read group
readgroup = illumina-pe-800
data = ${PE800}
autopairing
segment_placement = ---> <---
technology = solexa
EOF

        elif [ "$5" = "false" ]; then
            cat <<EOF >> $OUTPUT_FILE

#pe170 read group
readgroup = illumina-pe-170
data = ${PE170}
template_size = ${PE170_INSZ_MIN_MAX}
segment_placement = ---> <---
technology = solexa

#pe500 read group
readgroup = illumina-pe-500
data = ${PE500}
template_size = ${PE500_INSZ_MIN_MAX}
segment_placement = ---> <---
technology = solexa

#pe800 read group
readgroup = illumina-pe-800
data = ${PE800}
template_size = ${PE800_INSZ_MIN_MAX}
segment_placement = ---> <---
technology = solexa
EOF

        else
            printf "ERROR: Please set a valid value for the 5th argument: [true|false].\n\n"
            display_usage
            exit 1
        fi


    else
        printf "ERROR: File regarding the list of illumina pairs is not valid. Change the '$LIST_ILLUM_PE_MP_PATH' variable in
        the script to a valid file.\n"
        display_usage
        exit 1
    fi
fi









#############Mate pairs
if [ "$3" = "true" ]; then
    #check if is file
    if [ -f ${LIST_ILLUM_PE_MP_PATH} ]; then

        #loop to separate MP libraries by insert size to create the read group
        while read line
        do
            if [[ "$line" == *"MP2000"* ]] ; then
                MP2000="${MP2000}fastq::$line "

            elif [[ "$line" == *"MP5000"* ]] ; then
                MP5000="${MP5000}fastq::$line "

            else
                continue
            fi

        done < ${LIST_ILLUM_PE_MP_PATH}

        #check the auto estimate parameter
        if [ "$5" = "true" ] ; then
            cat <<EOF >> $OUTPUT_FILE
            
#mp2000 read group
readgroup = illumina-mp-2000
data = ${MP2000}
autopairing
segment_placement = <--- --->
technology = solexa

#mp5000 read group
readgroup = illumina-pe-5000
data = ${MP5000}
autopairing
segment_placement = <--- --->
technology = solexa
EOF

        elif [ "$5" = "false" ]; then
            cat <<EOF >> $OUTPUT_FILE
            
#mp2000 read group
readgroup = illumina-mp-2000
data = ${MP2000}
template_size = ${MP2000_INSZ_MIN_MAX}
segment_placement = <--- --->
technology = solexa

#mp5000 read group
readgroup = illumina-pe-5000
data = ${MP5000}
template_size = ${MP5000_INSZ_MIN_MAX}
segment_placement = <--- --->
technology = solexa
EOF

        else
            printf "ERROR: Please set a valid value for the 5th argument: [true|false].\n\n"
            display_usage
            exit 1
        fi

    else
        printf "ERROR: String regarding the list of illumina pairs is not valid. Change the '$LIST_ILLUM_PE_MP_PATH' variable in the script to a valid file.\n\n"
        display_usage
        exit 1
    fi
fi



###454 GENOMIC ########
if [ "$4" = "true" ]; then

        #check if path to 454 files exists
        if [ -d "$GENOMIC_454_PATH" ] ; then
            cat <<EOF >> $OUTPUT_FILE

#454 read group
readgroup = 454-genomic-data
data = ${GENOMIC_454_PATH}
technology = 454
EOF
        else
            printf "ERROR: String regarding the path for the 454 fastq files is not valid. Change the '$GENOMIC_454_PATH' variable in the script to a valid path.\n\n"
            display_usage
            exit 1
        fi
fi


##### ION TORRENT ###########
if [ "$6" = "true" ] && [ -d "$ION_TORRENT_PATH" ]; then

cat <<EOF >> $OUTPUT_FILE

#Ion torrent read group
readgroup = ion-torrent-data
data = ${ION_TORRENT_PATH}
technology = iontor
EOF
else
    printf "ERROR: String regarding the path for the IonTorrent fastq files is not valid. Change the '$ION_TORRENT_PATH' variable in the script to a valid path.\n\n"
    display_usage
    exit 1
fi




}




function settings(){
####PARAMETERS###

cat <<EOF >> $OUTPUT_FILE

#PARAMETERS
parameters = COMMON_SETTINGS -GE:not=$1:amm=no:mps=$2:kpmf=$3 -NW:cmrnl=warn \

EOF

#add parameters to force memory reduction
if [ -z "$4" ] || [ "$4" = "false" ]; then
    printf "Not forcing Mira to use less memory.\n"
cat <<EOF >> $OUTPUT_FILE
-SK:not=$1
EOF

else
cat <<EOF >> $OUTPUT_FILE
-SK:not=$1:mhpr=500:mhim=10000000
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
job = genome,denovo,accurate
EOF


####READ GROUPS
if [ -z "${10}" ]; then
    read_groups $1 $2 $3 $4 $5
elif [ "${10}" = "true" ] || [ [ "${10}" = "false" ]; then
    read_groups $1 $2 $3 $4 $5 "${10}"
else
    printf "\nERROR: Please provide a rigth value for the 10th parameter.\n\n"
    display_usage
    exit 1
fi

###PARAMETERS
if [ -z "$9" ]; then
    settings $6 $7 $8
elif [ "$9" = "true" ] || [ "$9" = "false" ]; then
    settings $6 $7 $8 $9
else
    printf "\nERROR: Please provide a rigth value for the 9th parameter.\n\n"
    display_usage
    exit 1
fi

