#!/bin/bash
if [ -z $1 ]; then
	printf "Error:Please provide a file that lists the path for the mapping files to process.\n\n"
	exit 1
fi

REFERENCE_ANNOTATION="/mnt/msa/otherProjects/populus-af/diffEpression/cufflinks-suite/Potrx01-genome.gtf"
REFERENCE_GENOME="/mnt/msa/otherProjects/populus-af/diffEpression/cufflinks-suite/filtered-plus500.fasta"
STRINGTIE="/opt/tools/stringtie-1.1.2/stringtie"
printf "####STEP 1####\n"
printf "Running individual StringTie assemblies for all the samples with the aid of the reference annotation...\n"
while read line 
do
	FILENAME=$line
	BASENAME=$(basename ${FILENAME})
	printf "Started sample $BASENAME...\n"
	LOG_FILE="logStringTie_${BASENAME}.txt"
#	mkdir assembly-${BASENAME}
	mkdir assembly2ndRound-${BASENAME}
#	$STRINGTIE $FILENAME -G $REFERENCE_ANNOTATION -o assembly-${BASENAME}/transcripts.gtf -A $PWD/assembly-${BASENAME}/gene-abundances.txt -C $PWD/assembly-${BASENAME}/cov-refs.gtf -p 10  &> $LOG_FILE
	#2nd round of stringtie [with -B and -e]
	$STRINGTIE $FILENAME -G mergedAssembly/merged.gtf -o $PWD/assembly2ndRound-${BASENAME}/transcripts.gtf -A $PWD/assembly2ndRound-${BASENAME}/gene-abundances.txt -C $PWD/assembly2ndRound-${BASENAME}/cov-refs.gtf -p 10 -B -e  &>> $LOG_FILE
	printf "Stringtie assembly for $BASENAME sample finished!!\n"
done < $1 

printf "Finished all stringtie assemblies!\n\n#####STEP 2####\nRunning merge of all individual assemblies with Cuffmerge..\n"


#for assembly in assembly*; do 
#	echo $assembly/transcripts.gtf >> individualAssemblies.txt; 
#done

#cuffmerge -g $REFERENCE_ANNOTATION -s $REFERENCE_GENOME -p 10 -o mergedAssembly individualAssemblies.txt &> log-cuffmerge.txt
#cuffmerge -g $REFERENCE_ANNOTATION  -p 10 -o mergedAssembly individualAssemblies.txt &> log-cuffmerge.txt
printf "Done.\n"
printf "Everything is finished."
