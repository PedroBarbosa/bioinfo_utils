#!/bin/bash
usage="$(basename "$0") combinationsFile[-h]. Script to get the common differential expressed genes across different methodologies and different combinations.

where:
    combinationsFile	Text file where each line represents a pair to be processed.
   
    Optional arguments:
    -h  show this help text"

if [ -z $1 ] ; then
	printf "Error. Few arguments.\n"
	echo "$usage"
	exit
fi

INTERCEPTION_SCRIPT="/mnt/msa/git_genosuber/genosuber/python-scripts/annotation/geneInterception.py"
RESULTS_FILE="interception-output.txt"
if [ -f "$RESULTS_FILE" ] ; then
	rm $RESULTS_FILE
fi

while read line 
do
	pair1=$(echo $line | awk '{print $1;}')	
	pair2=$(echo $line | awk '{print $2;}')
	printf "################################\n############$pair1 vs $pair2#############\n################################\n" > >(tee -a $PWD/$RESULTS_FILE)	
	printf "Processing DOWN regulated interceptions for $pair1 and $pair2 pair...\n" > >(tee -a $PWD/$RESULTS_FILE)
	ARRAY_DOWN=($(ls | grep "$pair1" | grep "$pair2" | grep -i "down"))
	printf "Number of files to process\t%i\n" "${#ARRAY_DOWN[@]}" > >(tee -a $PWD/$RESULTS_FILE)
	for file in "${ARRAY_DOWN[@]}"; do
		if [[ $file == *"edgeR"* ]] || [[ $file == *"edger"* ]]; then
			#Specific for the ouptu of my R script
			cat $file | cut -f1 -d "," | sed 's/\"//'g | tail -n +6 > "$file-ready"
 			
		elif [[ $file == *"cuffdiff"* ]] || [[ $file == *"cuffDiff"* ]]; then
			#Specific for the filterCuffDiff script output 
			cat $file | grep -v "^#" | cut -f1 -d " " > "$file-ready"
		fi
	done	
	
	FILES=($(ls | grep "ready"))
	PYTHON_RUN="python $INTERCEPTION_SCRIPT -g"
	for i in "${FILES[@]}"; do
		PYTHON_RUN="$PYTHON_RUN $i"
	done
#	echo $PYTHON_RUN 	
	$PYTHON_RUN &>> tmpFile.txt
	cat tmpFile.txt | grep -v "^(" | sed '/#/Q' &>> $RESULTS_FILE
	rm *ready
	rm unique*
	rm tmpFile.txt

	printf "Processing UP regulated interceptions for $pair1 and $pair2 pair...\n" > >(tee -a $PWD/$RESULTS_FILE)
        ARRAY_UP=($(ls | grep "$pair1" | grep "$pair2" | grep -i "up"))
        printf "Number of files to process\t%i\n\n\n" "${#ARRAY_UP[@]}"  > >(tee -a $PWD/$RESULTS_FILE)
        for file in "${ARRAY_UP[@]}"; do
                if [[ $file == *"edgeR"* ]] || [[ $file == *"edger"* ]]; then
                        #Specific for the ouptu of my R script
                        cat $file | cut -f1 -d "," | sed 's/\"//'g | tail -n +6 > "$file-ready"

                elif [[ $file == *"cuffdiff"* ]] || [[ $file == *"cuffDiff"* ]]; then
                        #Specific for the filterCuffDiff script output
                        cat $file | grep -v "^#" | cut -f1 -d " " > "$file-ready"
                fi
        done

	FILES=($(ls | grep "ready"))
	PYTHON_RUN="python $INTERCEPTION_SCRIPT -g"
        for i in "${FILES[@]}"; do
                PYTHON_RUN="$PYTHON_RUN $i"
	done
#       echo $PYTHON_RUN
        $PYTHON_RUN &>> tmpFile.txt
        cat tmpFile.txt | grep -v "^(" | sed '/#/Q' &>> $RESULTS_FILE
        rm *ready
        rm unique*
	rm tmpFile.txt

done < $1
