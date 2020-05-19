#!/bin/bash
if [ -z "$1" ] || [ -z "$2" ]; then
	printf "Error: Missing arguments!\n
	Usage:
		1st argument	File with the combination to filter. Each combination perl line, each sample separated by tab.
		2nd argument	Output file from cuffdiff.
		3rd argument	Optional. Value for the pvalue and qvalue to further filter the output file.\n"
		
	exit 1
fi
if [ ! -z "$3" ]; then
	pvalue=$3
fi

cuffdiff_file=$2
tempFile="tempFile.txt"
tempFile2="tempFile2.txt"
tempFile3="tempFile3.txt"
comparingFloatsScript="/mnt/msa/git_genosuber/genosuber/python-scripts/utilities/compare2floats.py"
while read line
do	
	
	sample1=$(echo $line | awk '{print $1;}')
	sample2=$(echo $line | awk '{print $2;}')
	outfileUp="${sample1}-${sample2}_cuffDiff_UP.txt"
	outfileDown="${sample1}-${sample2}_cuffDiff_DOWN.txt"
	if [ -f "$outfileUp" ]; then
		rm $outfileUp
	fi
	if [ -f "$outfileDown" ]; then
		rm $outfileDown
	fi
	header=$(head -1 $cuffdiff_file | cut -f1,5,6,8,9,10,12,13,14)
	printf "Processing $sample1 - $sample2 combination..\n"
	printf "Total number of differential expressed features:\t"
	cat $cuffdiff_file | grep -w "$sample1" | grep -w "$sample2" | grep -c "yes" 
#	myArray=($(cat $cuffdiff_file | cut -f1,5,6,8,9,10,12,13,14 | grep -w "$sample1" | grep -w "$sample2" | grep  "yes" | ))
	cat $cuffdiff_file | cut -f1,5,6,8,9,10,12,13,14 | grep -w "$sample1" | grep -w "$sample2" | grep  "yes" > $tempFile
	readarray myArray  < $tempFile
	up_genes=0
	down_genes=0
#	$(#myArray[@]})
#	printf "%s" "${myArray[@]}" > FODASS.txt 
	for line in "${myArray[@]}"
	do
		declare -a cols=($line)

#		COMPARING 2 FLOATS IN BASH, VERY SLOW
#		python $comparingFloatsScript $(echo ${cols[5]}) 0
#		if [ $(echo ${cols[5]} '>' 0 | bc -l) -eq 1 ]; then
#			let "up_genes += 1"	
#			echo $line >> $tempFile2	
#		elif [ $(echo ${cols[5]} '<' 0 | bc -l) -eq 1 ]; then
#			let "down_genes += 1"
#			echo $line >> $tempFile3

		if [ ${cols[5]} == "inf" ]; then #when fold change is inf
			let "up_genes += 1"
			echo $line >> $tempFile2
		
		elif [ ${cols[5]} == "-inf" ]; then
			let "down_genes +=1"
			echo $line >> $tempFile3

		else
			result=$(python $comparingFloatsScript $(echo ${cols[5]}) 0)
			if [ $result == "True" ]; then
				let "up_genes += 1"
	         		echo $line >> $tempFile2
			elif [ $result == "False" ] ; then
				let "down_genes +=1"
				echo $line >> $tempFile3
			else
				printf "LogFold change equal to 0. Something weird happened. Please check the script manually!\n"
				exit 1
			fi
		fi
	done

	printf "#Number of upregulated features of $sample2 over $sample1:\t$up_genes\n" > >(tee $PWD/$outfileUp)
	printf "#Number of downregulated features of $sample2 over $sample1:\t$down_genes\n" > >(tee $PWD/$outfileDown)
	
	#insert header
	echo "#$header" | cut -f1,4,5,6,7,8 >> $outfileUp	
	echo "#$header" | cut -f1,4,5,6,7,8 >> $outfileDown
	#sort file by pvalue and write the most significant 50 genes to each file
	sort -k7 -g $tempFile2 | cut -f1,4,5,6,7,8 -d " " | head -50 >> $outfileUp
	sort -k7 -g $tempFile3 | cut -f1,4,5,6,7,8 -d " " | head -50 >> $outfileDown
	
	#all diffExp genes for intercpetion analysis
	sort -k7 -g $tempFile2 | cut -f1,4,5,6,7,8 -d " " > all-${outfileUp}		
	sort -k7 -g $tempFile3 | cut -f1,4,5,6,7,8 -d " " > all-${outfileDown}
	#remove temp files
	rm $tempFile
	rm $tempFile2
	rm $tempFile3
done < $1
