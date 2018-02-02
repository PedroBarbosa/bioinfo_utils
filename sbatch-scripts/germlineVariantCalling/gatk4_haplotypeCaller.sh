#!/bin/bash
display_usage(){
 printf "Script to run Haplotype caller pipeline in lobo for a bunch of bam files. Be aware that this script already implements the latest innovations regarding
   the multisample variant calling in GVCF mode and ImportGenomicsDB.\n
 Usage:.
    -1st argument must be the file containing the alignment files. A '.bai' index for each bam file must exist in each bam directory.
    -2nd argument must be the name of the final ouptut directory.
    -3rd argument must refer to the type of analysis in hand. Values:[WGS|targeted].
    -4th argument must be the fasta reference for which the reads were aligned. If '-' is set, the existing hg38 version in lobo will be used. Be aware that a fasta index file (.fai) must also be present in the fasta directory.
    -5th argument is optional. If set to false, GNU parallel will be disabled to run the set of samples provided in the 1st argument. Options: [true|false]. Default: true, GNU parallel is used to parallelize the job.
    -6th argument must be provided when the data comes from a targeted experiment. Refers to the target intervals in bed format. Use '-' to skip this argument.
    -7th argument is optional. If set, refers to a known variants file (e.g dbsnp) with IDs.  Its purpose is to annotate our variants with the corresponding reference ID.\n"
}

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ]; then
    printf "Error. Please set the required parameters of the script\n"
    display_usage
    exit 1
fi

####CHECK BAM INPUT####
if [ ! -f "$1" ]; then
    echo "File "$1" not found!"
    exit 1
fi
while read line
do
    if [ ! -f "$line" ]; then
         printf "Error. $line doesn't exist. Please check more carefully files passed in the 1st argument.\n"
         display_usage
         exit 1
    fi
done < "$1"
BAM_DATA=$(readlink -f "$1")

####SCRATCH WORKDIR####
OUTDIR=$(readlink -f "$2")
WORKDIR="/home/pedro.barbosa/scratch/gatk"
if [ ! -d $WORKDIR ];then
    mkdir $WORKDIR
fi
if [ ! -d $OUTDIR ];then
    mkdir $OUTDIR
fi

##REFERENCE##
if [ "$4" = "-" ]; then
    REF="/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/WholeGenomeFasta/genome.fa"
elif [ ! -f "$4" ]; then
    printf "Please provide a valid fasta file in th 4th argument.\n"
    display_usage
    exit 1
elif [ ! -f "${4}.fai" ]; then
    printf "Fasta index ${4}.fai not found in the reference directory. Please create one with samtools faidx.\n"
    display_usage
    exit
else
    REF="$4"
fi

###PARALLEL####
if [ -z "$5" ] || [ "$5" = "true" ]; then
    NODES=1 #2 #1
    NTASKS=5 #10 #4
    CPUS=5 #8 #10
    JAVA_Xmx="--java-options '-Xmx45G'"
    PARALLEL=true
elif [ "$5" = "false" ]; then
    NODES=1
    NTASKS=1
    CPUS=40
    JAVA_Xmx="--java-options '-Xmx240G'"
    PARALLEL=false
else
    printf "Please set a valid value for the 5th argument\n"
    display_usage
    exit 1
fi
##CMD##
CMD="gatk ${JAVA_Xmx} HaplotypeCaller --reference=$REF -ERC GVCF"

###MODE####
analysis=("WGS" "targeted")
if [[ !  " ${analysis[@]} " =~ " ${3} " ]]; then
    printf "Please set a valid value for the 3th argument.\n"
    display_usage
    exit 1
elif [ "$3" != "WGS" ]; then
    if [ -z "$6" ]; then
        printf "When exome or targeted sequencing, you should provide the target regions in the 6th argument as a single bed file.\n"
        display_usage
        exit 1
    elif [[ ${6: -4} == ".bed"  ]]; then
        BED=$(readlink -f "$6")
        CMD="$CMD --intervals="$BED" --interval-padding 50"
    else
        printf "Invalid bed input for target regions.\n"
        display_usage
        exit 1
    fi  
fi


##Dbsnp file###

cat > $WORKDIR/haplotypeCaller_GVCF.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=gatk_hc
#SBATCH --time=72:00:00
#SBATCH --mem=240G
#SBATCH --nodes=$NODES
#SBATCH --ntasks=$NTASKS
#SBATCH --cpus-per-task=$CPUS
#SBATCH --image=broadinstitute/gatk:latest
#SBATCH --workdir=$WORKDIR
#SBATCH --output=$WORKDIR/%j_gatk_hc.log

timestamp() {
    date +"%Y-%m-%d  %T"
}
export -f timestamp
scratch_out=$WORKDIR/\$SLURM_JOB_ID
mkdir \$scratch_out
cd \$scratch_out
srun="srun -N1 -n1 --slurmd-debug 3"
parallel="parallel --delay 0.2 -j \$SLURM_NTASKS  --env timestamp --joblog parallel.log --resume-failed"
echo "\$(timestamp) -> Job started!"
if [ $PARALLEL = "true" ]; then
    cat $BAM_DATA | \$parallel '\$srun shifter \$CMD -I={} -O={/.}.g.vcf'
else
    unique_samples=()
    for j in \$(find $BAM_DATA -exec cat {} \; );do
        printf "\$(timestamp): Processing \$(basename \$j) file.\n"
        i=\$(basename \$j)
        out=\$(echo \$i | cut -f1 -d "_")
        if [[ " \${unique_samples[@]} " =~ " \${out} " ]]; then
            printf "\$(timestamp): Warning, duplicate sample ID (\${out}) found after splitting by first '_'. Will use the filename instead without the extension."
            out="\${i%.*}"
        fi
        unique_samples+=(\${out})
        \$srun shifter \$CMD -I=\$j -O=\${out}.g.vcf
        echo "\${out}\t$\j" >> aux_sampleName.map
        printf "\$(timestamp): Done!\n"
    done
fi
echo "\$(timestamp) -> Haplotype caller in GVCF for all samples has finished."
echo "\$(timestamp) -> Creating GenomicsDB datastore to combine multiple gvcf data."
##TRY ALSO FOR ONE SINGLE INTERVALS FILE###
CMD_GenomicsDB="./gatk GenomicsDBImport --sample-name-map aux_sampleName.map --batch-size 5 -R $REF --genomicsdb-workspace-path workspace"
if [ "$3" = "WGS" ]; then
    for i in $(seq 1 22); do echo -e "\$(timestamp) Consolidating gVCF for chromosome \${i}"; \$srun shifter \$CMD_GenomicsDB -L chr\${i}; done
    echo -e "\$(timestamp) Consolidating gVCF for sexual chromosomes and mithocondria."
    \$srun shifter \$CMD_GenomicsDB -L chrM
    \$srun shifter \$CMD_GenomicsDB -L chrY
    \$srun shifter \$CMD_GenomicsDB -L chrX
else
    for j in \$(find $BED -exec cat {} \; );do
        echo -e "\$(timestamp) Consolidating gVCF for \$(echo \$j | cut -f1,2,3) interval."
        \$srun shifter \$CMD_GenomicsDB -L \$(echo \$j | cut -f1):\$((\$(echo \$j | cut -f2) + 1))-\$(echo \$j | cut -f3)
    done
echo "\$(timestamp) -> Done! Genotyping gVCFs stored in the GenomicsDB database."
\$srun shifter gatk GenotypeGVCFs --variant gendb://workspace -R $REF --output joint.vcf.gz
echo -e "\$(timestamp) -> Finished job."
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
echo -e "\$(timestamp) -> All done!"
EOL
sbatch $WORKDIR/haplotypeCaller_GVCF.sbatch
sleep 1 
cd $WORKDIR
mv haplotypeCaller_GVCF.sbatch $(ls -td -- */ | head -n 1)
