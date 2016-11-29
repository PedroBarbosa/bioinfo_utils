#/bin/bash
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
	printf "Few arguments\n"
	printf "Usage:

		-1st argument	File listing individual assemblies to process[gtf file produced by stringTie]
		-2nd argument	Merged transcriptome file[gtf file produced by cuffmerge with all samples]
		-3rd argument	File listing gene abundances produced by each individual StringTie assembly\n\n"
	exit 1
fi


total_genes=$(cat $2 | grep -v "#" | cut -f9 | cut -f1 -d ";" | sort | uniq | wc -l)
total_isoforms=$(cat $2 | grep -v "#" | cut -f9 | cut -f2 -d ";" | sort | uniq | wc -l)
lines=1
while read line
do	
#	split_sample_path=(${line//\// }) #split by "/"
#	sample_name=${split_sample_path[-2]}
#	sample_name=$(basename $line | sed 's/\(_*\).*/\1/')
	sample_name=$line	

#	unique_features=$(tail -n +2 $line | cut -f1 | sort | uniq | wc -l)
#	all_features=$(tail -n +2 $line | cut -f1,8,9 | sort | uniq |  wc -l)
#	expressed_features=$(tail -n +2 $line | cut -f1,8 | sort | uniq | cut -f2 | grep -cwv "0")  
#	duplicate_features=$(tail -n +2 $line | cut -f1,10 | cut -f1 |sort | uniq -cd | wc -l)

	#GENE ANALYSIS#
	unique_genes_stringtie=$(cat $line | grep -w "transcript" | cut -f9 | cut -f1 -d ";" | sort | uniq | wc -l)
	#unique_genes_reference=$(cat $line | grep -w "transcript" | cut -f9 | cut -f4 -d ";" | sort | uniq | wc -l)
	#unique_genes_both=$(cat $line | grep -w "transcript" | cut -f9 | cut -f1,4 -d ";" | sort | uniq | wc -l)
	#unique_genes_fromReference=$(cat $line | grep -w "transcript" | cut -f9 | cut -f5 -d ";" | grep "ref_gene_name" | sort | uniq | wc -l)
	gene_abundance_file=$(sed -n "${lines}p" < $3)
	unique_genes_abundance_file=$(tail -n +2 $gene_abundance_file | cut -f1 | sort | uniq | wc -l)
        let "lines +=1"


	#ISOFORM ANALYSIS#
	unique_isoforms=$(cat $line | grep -v "#" | grep -w "transcript" | cut -f9 | cut -f2 -d ";" | sort | uniq | wc -l)
	#expressed_isoforms=$(cat $line | grep -v "#" | grep -w "transcript" | cut -f9 | cut -f2,6,7,8 -d ";" | grep -c "FPKM")
        expressed_isoforms=$(cat $line | grep -v "#" | grep -w "transcript" | cut -f9 | cut -f2,3,4,5 -d ";" | grep -c "FPKM")	
 
	printf "###SAMPLE $sample_name###\n"
	printf "Number of genes tracked:\t$total_genes\n"
	printf "Number of isoforms tracked:\t$total_isoforms\n\n"
	printf "Number of unique stringTie genes expressed in sample:\t$unique_genes_stringtie\n"
#	printf "Number of unique reference [cuffmerge IDs] genes expressed in sample:\t$unique_genes_reference\n"
#	printf "Number of unique reference and stringTie genes expressed in sample:\t$unique_genes_both\n"
#	printf "Number of genes expressed mapping reference genes from original GTF:\t$unique_genes_fromReference\n"
	printf "Number of unique genes expressed in the gene abundances file:\t$unique_genes_abundance_file\n\n"
#	printf "Number of unique reference and StringTie isoforms present in sample:\t$unique_isoforms\n"
	printf "Number of unique isoforms present in sample:\t$unique_isoforms\n"
#	printf "Number of unique reference and StringTie isoforms expressed in sample:\t$expressed_isoforms\n\n\n\n"
	printf "Number of unique StringTie isoforms expressed in sample:\t$expressed_isoforms\n\n\n\n"

#	printf "Duplicated features in file:\t%s\n\n\n" "$(tail -n +2 $line | cut -f1,10 | cut -f1 |sort | uniq -cd)"
#	total=$(($unique_features + $duplicate_features)) #)cat in.fasta | wc -l) / 4

#	printf "Number to check:\t$total\n\n\n"

done < $1
