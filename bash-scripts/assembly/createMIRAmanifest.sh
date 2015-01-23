#!/bin/bash

display_usage(){
 printf "Script to automatically generate the manifest file to input in Mira for de novo assembly projects. Paired end
 reads mus be in the traditional 'forward-reverse' orientation while mate pair must come in the 'reverse-forward' orientation'.\n\n
 Usage:
    -1st argument must be the project name to use in Mira.The name will be the prefix for the manifest file generated.
    -2nd argument must be a flag to use paired end reads to generate the command.Available option: [true|false].
    -3rd argument must be a flag to use mate pair reads to generate the command.Available option: [true|false].
    -4th argument must be a flag to use 454 reads reads to generate the command.Available option: [true|false].
    -5th argument must be a flag to Mira auto estimate insert sizes of libraries [true|false].
"
}


####################VARIABLES####################
#Project name
PROJECT_NAME="$1"

#Output file
OUTPUT_FILE="$PROJECT_NAME-manifest.txt"

#Path to the genomic 454 sequencing FASTA files - All files of the directory will be used so we can provide a path in Mira
GENOMIC_454_PATH="/mnt/msa/BIOCANT/genomic-data/SFF_genom/FASTA_FILES/"

#Path to the file that lists Illumina files wit the threshold Q20L80/Q20L40
LIST_ILLUM_PE_MP_PATH="/mnt/msa/workflow_scripts/LIST_FILES/listFiles_Q20L80-PE_Q20L20-MP.txt"

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

##############################START APPENDING GENERAL SETTINGS TO OUTPUT FILE##########
cat <<EOT >> $OUTPUT_FILE
#Name to the assembly and purpose of the job
project = $PROJECT_NAME
job = genome,denovo,accurate
EOT



###################DEFINING READ GROUPS#############################


#If assemble paired end
if [ "$2" = "true" ]; then
    #check if is file
    if [ -f ${LIST_ILLUM_PE_MP_PATH} ]; then

        #loop to separate PE libraries by insert size to create the read group
        while read line
        do
            if [[ "$line" == *"PE170"* && "$line" =~ *"Single"* ]] ; then
                PE170="$PE170 $line"

            elif [[ "$line" == *"PE500"* && "$line" =~ *"Single"* ]] ; then
                PE500="$PE500 $line"

            elif [[ "$line" == *"PE800"* && "$line" =~ *"Single"* ]] ; then
                PE800="$PE800 $line"

            else
                continue
            fi

        done < $LIST_ILLUM_PE_MP_PATH

        #check the auto estimate parameter
        if [ "$5" = "true" ] ; then
            cat <<EOT>> $OUTPUT_FILE
                #pe170 read group
                readgroup = illumina-pe-170
                data = ${PE170}
                autopairing
                segment_placement = ---> <---
                segment_naming = solexa

                #pe500 read group
                readgroup = illumina-pe-500
                data = ${PE500}
                autopairing
                segment_placement = ---> <---
                segment_naming = solexa

                #pe800 read group
                readgroup = illumina-pe-800
                data = ${PE800}
                autopairing
                segment_placement = ---> <---
                segment_naming = solexa
            EOT

        elif [ "$5" = "false" ]; then
            cat <<EOT>> $OUTPUT_FILE
                #pe170 read group
                readgroup = illumina-pe-170
                data = ${PE170}
                template_size ${PE170_INSZ_MIN_MAX}
                segment_placement = ---> <---
                segment_naming = solexa

                #pe500 read group
                readgroup = illumina-pe-500
                data = ${PE500}
                template_size ${PE500_INSZ_MIN_MAX}
                segment_placement = ---> <---
                segment_naming = solexa

                #pe800 read group
                readgroup = illumina-pe-800
                data = ${PE800}
                template_size ${PE800_INSZ_MIN_MAX}
                segment_placement = ---> <---
                segment_naming = solexa
            EOT
        else
            printf "Please set a valid value for the 5th argument: [true|false].\n\n"
            display_usage
            exit 1
        fi
            data =

    else
        printf "File regarding the list of illumina pairs is not valid. Change the '$LIST_ILLUM_PE_MP_PATH' variable in
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
            if [[ "$line" == *"MP2000"* && "$line" =~ *"Single"* ]] ; then
                MP2000="$MP2000 $line"

            elif [[ "$line" == *"MP5000"* && "$line" =~ *"Single"* ]] ; then
                MP5000="$MP5000 $line"

            else
                continue
            fi

        done < $LIST_ILLUM_PE_MP_PATH

        #check the auto estimate parameter
        if [ "$5" = "true" ] ; then
            cat <<EOT>> $OUTPUT_FILE
                #mp2000 read group
                readgroup = illumina-mp-2000
                data = ${MP2000}
                autopairing
                segment_placement = <--- --->
                segment_naming = solexa

                #mp5000 read group
                readgroup = illumina-pe-5000
                data = ${MP5000}
                autopairing
                segment_placement = <--- --->
                segment_naming = solexa

            EOT

        elif [ "$5" = "false" ]; then
            cat <<EOT>> $OUTPUT_FILE
                #mp2000 read group
                readgroup = illumina-mp-2000
                data = ${MP2000}
                template_size ${MP2000_INSZ_MIN_MAX}
                segment_placement = <--- --->
                segment_naming = solexa

                #mp5000 read group
                readgroup = illumina-pe-5000
                data = ${MP5000}
                template_size ${MP2000_INSZ_MIN_MAX}
                segment_placement = <--- --->
                segment_naming = solexa

            EOT
        else
            printf "Please set a valid value for the 5th argument: [true|false].\n\n"
            display_usage
            exit 1
        fi

    else
        printf "File regarding the list of illumina pairs is not valid. Change the '$LIST_ILLUM_PE_MP_PATH' variable in
        the script to a valid file.\n"
        display_usage
        exit 1
    fi
fi



###454 GENOMIC ########
if [ "$4" = "true" ]; then

        #check if path to 454 files exists
        if [ -d "$GENOMIC_454_PATH" ] ; then
            cat <<EOT>> $OUTPUT_FILE
                #454 read group
                readgroup = 454-genomic-data
                data = FASTA_FILES/
                technology = 454
            EOT
        else
            printf "File regarding the path for the 454 fasta files is not valid. Change the '$GENOMIC_454_PATH' variable in
        the script to a valid path.\n"
            display_usage
            exit 1
        fi
fi




####PARAMETERS###
#parameters = COMMON SETTINGS -GE:not=
