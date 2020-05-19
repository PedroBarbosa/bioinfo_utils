#!/bin/bash
#This is a script that creates a tabular file with the existing files present in a directory representing the output of the 'estimate insert size' script present in the SSPACE software. It simplifies the visualization of the insert size estimation for several libraries
DIR=./MP
printf "File\tFF_pairs\tFF_ins\tFF_stdv\t\tRR_pairs\tRR_ins\tRR_stdv\t\tRF_pairs\tRF_ins\tRF_stdv\t\tFR_pairs\tFR_ins\tFR_stdv\n" > MP_insert_estimation.txt

FF=false
RR=false
RF=false
FR=false
write2file=false
previous_file=""
for file in $DIR/*
do	
	if [ "$write2file" = "true" ] ; then
			printf "$previous_file\t$pairs_FF\t$insert_FF\t$stdev_FF\t\t$pairs_RR\t$insert_RR\t$stdev_RR\t\t$pairs_RF\t$insert_RF\t$stdev_RF\t\t$pairs_FR\t$insert_FR\t$stdev_FR\n" >> MP_insert_estimation.txt
	fi
	write2file=true
	echo $file
	previous_file=$file
	while read line
	do		
		
		##FF###
		if [[ "$line" == *"FF"* ]]; then
			FF=true
			RR=false
			RF=false
			FR=false	
		elif [[ "$line" == *"RR"* ]]; then
			RR=true
			FF=false
			RF=false
			FR=false		
		elif [[ "$line" == *"RF"* ]]; then
			RF=true
			FF=false
			RR=false
			FR=false
		elif [[ "$line" == *"FR"* ]]; then
			FR=true
			FF=false
			RR=false
			RF=false
		fi
		
		
		
		if [[ "$line" == *"pairs"* && "$FF" = "true" ]]; then
			pairs_FF=$(echo $line | egrep -o '[0-9]+')
		
		elif [[ "$line" == *"median insert"* && "$FF" = "true" ]]; then
			insert_FF=$(echo $line | egrep -o '[0-9]+')
		
		elif [[ "$line" == *"stdev"* && "$FF" = "true" ]]; then
			stdev_FF=$(echo $line | egrep -o '[0-9.0-9]+')
		fi
	
		##RR##
		if [[ "$RR" = "true" && "$line" == *"pairs"* ]]; then
			pairs_RR=$(echo $line | egrep -o '[0-9]+')
		
		elif [[ "$RR" = "true" && "$line" == *"median insert"* ]]; then
			insert_RR=$(echo $line | egrep -o '[0-9]+')
		
		elif [[ "$RR" = "true" && "$line" == *"stdev"* ]]; then
			stdev_RR=$(echo $line | egrep -o '[0-9.0-9]+')
		fi
		
		
		
		##RF###
		if [[ "$RF" = "true" && "$line" == *"pairs"* ]]; then
			pairs_RF=$(echo $line | egrep -o '[0-9]+')
		
		elif [[ "$RF" = "true" && "$line" == *"median insert"* ]]; then
			insert_RF=$(echo $line | egrep -o '[0-9]+')
		
		elif [[ "$RF" = "true" && "$line" == *"stdev"* ]]; then
			stdev_RF=$(echo $line | egrep -o '[0-9.0-9]+')
		fi
		
		
		##FR###
		if [[ "$FR" = "true" && "$line" == *"pairs"* ]]; then
			pairs_FR=$(echo $line | egrep -o '[0-9]+')
		
		elif [[ "$FR" = "true" && "$line" == *"median insert"* ]]; then
			insert_FR=$(echo $line | egrep -o '[0-9]+')
		
		elif [[ "$FR" = "true" && "$line" == *"stdev"* ]]; then
			stdev_FR=$(echo $line | egrep -o '[0-9.0-9]+')
		fi


	done < $file
done


