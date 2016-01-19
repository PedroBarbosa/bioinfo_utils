#/bin/bash
if [ -z "$1" ] ; then
	printf "Please add the a file with the list of files to process\n"
	printf "This script requires gene or isoforms tracking files produced by cufflinks, not the GTF file\n"
	exit 1
fi


while read line
do	
#	declare -A feature_fpkm
	split_sample_path=(${line//\// }) #split by "/"
	sample_name=${split_sample_path[-2]}
	

	unique_features=$(tail -n +2 $line | cut -f1 | sort | uniq | wc -l)
	all_features=$(tail -n +2 $line | cut -f1,10 | sort | uniq |  wc -l)
	expressed_features=$(tail -n +2 $line | cut -f1,10 | sort | uniq | cut -f2 | grep -cwv "0")  
	duplicate_features=$(tail -n +2 $line | cut -f1,10 | cut -f1 |sort | uniq -cd | wc -l)
	featuresExpressedFromOriginalGTF=$(tail -n +2 $line | cut -f1,10 | grep -v "CUFF" | sort | uniq | cut -f2 | grep -cwv "0")
#tail -n +2 $line | cut -f1,10 | cut -f1 |sort | uniq -cd
	printf "###SAMPLE $sample_name###\n"
	printf "Number of features tracked:\t$all_features\n"
	printf "Number expressed features:\t$expressed_features\n"
	printf "Number of features expressed mapping reference genes from original GTF:\t$featuresExpressedFromOriginalGTF\n"
	printf "Number of unique features in file:\t$unique_features\n"
	printf "Number of duplicate features in file:\t$duplicate_features\n"
	printf "Duplicated features in file:\t%s\n\n\n" "$(tail -n +2 $line | cut -f1,10 | cut -f1 |sort | uniq -cd)"
#	total=$(($unique_features + $duplicate_features)) #)cat in.fasta | wc -l) / 4
#	printf "Number to check:\t$total\n\n\n"


#	while read line
#	do

#		feature_fpkm[$(echo $line | cut -f1 -d " ")]=$(echo $line | cut -f10 -d " ")
#	done < $line

#	for feature in "${!feature_fpkm[@]}"; do
#		if [ ${feature_fpkm[$feature]} != "0" ]; then
#			let "expressed_genes+=1"		
#		fi
#		let "total_genes+=1"
#	done
	
#	printf "Number of genes tracked in $sample_name:\t$total_genes\n"
#	printf "Number of genes tracked with expression observed in $sample_name:\t$expressed_genes\n\n"

done < $1
