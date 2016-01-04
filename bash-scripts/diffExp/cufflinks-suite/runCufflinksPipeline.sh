#!/bin/bash
if [ -z $1 ]; then
        printf "Error:Please provide a file that lists the path for the mapping files to process.\n\n"
        exit 1
fi

REFERENCE_ANNOTATION="/mnt/msa/otherProjects/populus-af/diffEpression/cufflinks-suite/Potrx01-genome.gtf"
REFERENCE_GENOME="/mnt/msa/otherProjects/populus-af/diffEpression/cufflinks-suite/filtered-plus500.fasta"

printf "####STEP 1####\n"
printf "Running individual Cufflinks assemblies for all the samples without the aid of the reference annotation...\n"
while read line
do
        FILENAME=$line
        BASENAME=$(basename ${FILENAME})
        printf "Started sample $BASENAME...\n"
        LOG_FILE="logCufflinks_${BASENAME}.txt"
#       cufflinks -o assembly-${BASENAME} -p 10 -u --library-type fr-firststrand $FILENAME &> $LOG_FILE
        cufflinks -g $REFERENCE_ANNOTATION -b $REFERENCE_GENOME -o assembly-${BASENAME} -p 10 -u --library-type fr-firststrand $FILENAME &> $LOG_FILE
        printf "Cufflinks assembly for $BASENAME sample finished!!\n"
done < $1

printf "Finished all Cufflinks assemblies!\n\n#####STEP 2####\nRunning merge of all individual assemblies with Cuffmerge..\n"


for assembly in assembly*; do
        echo $assembly/transcripts.gtf >> individualAssemblies.txt;
done

cuffmerge -g $REFERENCE_ANNOTATION -s $REFERENCE_GENOME -p 10 -o mergedAssembly individualAssemblies.txt &> log-cuffmerge.txt
#cuffmerge -g $REFERENCE_ANNOTATION  -p 10 -o mergedAssembly individualAssemblies.txt &> log-cuffmerge.txt
printf "Cuffmerged finished!!\n\n#####STEP 3####\nRunning CuffDiff to find differential expression between conditions..\n"
printf "Running cuffDiff all vs all..\n"
cuffdiff -o all-vs-all -b $REFERENCE_GENOME -p 10 -u --library-type fr-firststrand -L c1,c2,c3,c4,c5,c3c4c5 mergedAssembly/merged.gtf /mnt/msa/otherProjects/populus-af/mapping/star/stem-samples/C1_RNASeq_Pop_L1_UMR_noS.bam,/mnt/msa/otherProjects/populus-af/mapping/star/stem-samples/C1_RNASeq_Pop_L2_UMR_noS.bam /mnt/msa/otherProjects/populus-af/mapping/star/stem-samples/C2_RNASeq_Pop_L1_UMR_noS.bam,/mnt/msa/otherProjects/populus-af/mapping/star/stem-samples/C2_RNASeq_Pop_L2_UMR_noS.bam /mnt/msa/otherProjects/populus-af/mapping/star/stem-samples/C3_RNASeq_Pop_L1_UMR_noS.bam,/mnt/msa/otherProjects/populus-af/mapping/star/stem-samples/C3_RNASeq_Pop_L2_UMR_noS.bam /mnt/msa/otherProjects/populus-af/mapping/star/stem-samples/C4_RNASeq_Pop_L1_UMR_noS.bam,/mnt/msa/otherProjects/populus-af/mapping/star/stem-samples/C4_RNASeq_Pop_L2_UMR_noS.bam /mnt/msa/otherProjects/populus-af/mapping/star/stem-samples/C5_RNASeq_Pop_L1_UMR_noS.bam,/mnt/msa/otherProjects/populus-af/mapping/star/stem-samples/C5_RNASeq_Pop_L2_UMR_noS.bam /mnt/msa/otherProjects/populus-af/mapping/star/stem-samples/C3_RNASeq_Pop_L1_UMR_noS.bam,/mnt/msa/otherProjects/populus-af/mapping/star/stem-samples/C3_RNASeq_Pop_L2_UMR_noS.bam,/mnt/msa/otherProjects/populus-af/mapping/star/stem-samples/C4_RNASeq_Pop_L1_UMR_noS.bam,/mnt/msa/otherProjects/populus-af/mapping/star/stem-samples/C4_RNASeq_Pop_L2_UMR_noS.bam,/mnt/msa/otherProjects/populus-af/mapping/star/stem-samples/C5_RNASeq_Pop_L1_UMR_noS.bam,/mnt/msa/otherProjects/populus-af/mapping/star/stem-samples/C5_RNASeq_Pop_L2_UMR_noS.bam &> log-cuffdiff.txt

printf "Done.\n"
printf "Everything is finished."

