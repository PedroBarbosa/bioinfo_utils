#!/bin/bash
#Script to run bwa for each library with default parameters. One sorted and indexed bam per library will be produced.
display_usage() {
echo '1st argument must be the file. In this file the Paired end libraries need to be first
2nd argument must be a flag true/false to use paired end reads to generate the command
3rd argument must be a flag true/false to use mate pair reads to generate the command
4th argument must be the prefix for the reference indexed database
5th argument must be the number of threads to use
6th argument is optional. Must be string to add to the read group parameter.Add the LB, PL and PU fields, respectively,comma separated. Ex: LB:lib,PL:instrument,PU:plataformUnit. The SM and ID fields will be automatically added based on the sample basename.'
}


exec="bwa mem"
#check if required arguments are there and display usage message
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ]; then
        printf "Please provide the arguments required for the script.\n\n"
        display_usage
        exit 1
fi

if [ -z "$6" ]; then
        readGroup="false"
else
        readGroup="true"
        read_group_general=$(echo "${6//,/\t}")
fi

index_database="$4"
threads="-t $5"
numb_samples=0
first_pair=true
second_pair=false
matepairFlag=false
base_command="$exec $threads $index_database "


while read line
do

        #remove leading and trailing spaces
        #filename=$(echo $line | sed -e 's/^ *//' -e 's/ *$//')
        filename=$line

        #check if reached the MATE PAIR samples, and if not supposed to map them, break the loop
        if [[ "$filename" == *"MP"* && "$3" = "true" ]];then
                matepairFlag=true

        elif [[ "$filename" == *"MP"* && "$3" = "false" ]];then
                break
        fi

        bam_basename=$(basename $filename | sed 's/\(.*\)_.*/\1/')
        #echo $bam_basename


        #PAIRED-END###
        if [ -f "$filename" -a "$first_pair" = "true" -a "$2" = "true" -a "$matepairFlag" = "false" ]; then
                pair1="$filename"
                first_pair=false
                second_pair=true
        elif [ -f "$filename" -a  "$second_pair" = "true" -a "$2" = "true" -a "$matepairFlag" = "false" ]; then

                pair2="$filename"
                #bam_basename=$(basename $pair2 | cut -d "_" -f1)
                first_pair=true
                second_pair=false
                let "numb_samples += 1"
                printf "\nRunning bwa mem aligner for library ${bam_basename} [sample ${numb_samples}] .."
                bam_file="${bam_basename}.bam"
                if [ "$readGroup" = "true" ];then
                        command="$base_command -R @RG\tID:${bam_basename}_id\tSM:${bam_basename}\t$read_group_general $pair1 $pair2"
                else
                        command="$base_command $pair1 $pair2"
                fi
                command_view="samtools view -Sbh -"
                command_sort="samtools sort ${bam_file}"
                command_index="samtools index"
		        printf "##CMD##:\n$command | $command_view > ${bam_file}\n"
      		    $command | $command_view > ${bam_file} 2>> ./stderr.txt
                printf "Done!! Sorting bam file..\n##CMD##\n$command_sort > ${bam_file/.bam/_sorted.bam}\n\n"
                $command_sort > ${bam_file/.bam/_sorted.bam} 2>> ./stderr.txt
                $command_index ${bam_file/.bam/_sorted.bam} 2>> ./stderr.txt
#                rm $bam_file

        elif [ ! -f "$filename" -a "$matepairFlag" = "false"  ]; then
                echo "$filename" is not a file
                exit 1
        fi


        ##MATE-PAIR$####
        if [ "$matepairFlag" = "true" -a -f "$filename" -a "$first_pair" = "true" -a "$3" = "true" ]; then
                pair1="$filename"
                first_pair=false
                second_pair=true

        elif [ "$matepairFlag" = "true" -a -f "$filename" -a  "$second_pair" = "true" -a "$3" = "true" ]; then


                pair2="$filename"
                #bam_basename=$(basename $pair2 | cut -d "_"
                first_pair=true
                second_pair=false
                let "numb_samples += 1"
                printf "\nRunning bwa mem aligner for library ${bam_basename} [sample ${numb_samples}] ..\n"

                bam_file="${bam_basename}.bam"
                if [ "$readGroup" = "true" ];then
                        command="$base_command -R @RG\tID:${bam_basename}_id\tSM:${bam_basename}\t$read_group_general $pair1 $pair2"
                else
                        command="$base_command $pair1 $pair2"
                fi

                command_view="samtools view -Sbh -"
                command_sort="samtools sort ${bam_file}"
                command_index="samtools index"
                printf "##CMD##:\n$command | $command_view > ${bam_file}\n"
                $command | $command_view > ${bam_file} 2>> ./stderr.txt
	            printf "Done!! Sorting bam file..\n##CMD##\n$command_sort > ${bam_file/.bam/_sorted.bam}\n\n"
                $command_sort > ${bam_file/.bam/_sorted.bam} 2>> ./stderr.txt
                $command_index ${bam_file/.bam/_sorted.bam} 2>> ./stderr.txt
#                rm $bam_file

        elif [ ! -f "$filename" -a "$matepairFlag" = "true" ]; then
                echo "$filename" is not a file
                continue
        fi

done <$1
