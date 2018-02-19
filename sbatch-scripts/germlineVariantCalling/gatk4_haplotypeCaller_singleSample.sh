#!/bin/bash
display_usage(){
 printf "Script to run Haplotype caller pipeline in single sample mode in lobo for a bunch of bam files. One VCF file per sample will be produced via HaplotypeCallerSpark utility.\n
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
    CPUS=7 #8 #10
    JAVA_Xmx="--java-options '-Xmx45G'"
    PARALLEL=true
elif [ "$5" = "false" ]; then
    NODES=1
    NTASKS=1
    CPUS=40
    JAVA_Xmx="--java-options '-Xmx240G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true'"
    PARALLEL=false
else
    printf "Please set a valid value for the 5th argument\n"
    display_usage
    exit 1
fi
##CMD##
ref_2bit="$(basename $REF)"
CMD="gatk ${JAVA_Xmx} HaplotypeCallerSpark --TMP_DIR=/home/pedro.barbosa/scratch --reference ${ref_2bit%.*}.2bit"
#CMD="gatk ${JAVA_Xmx} HaplotypeCaller --TMP_DIR=/home/pedro.barbosa/scratch --reference $REF -ERC GVCF"
#createOutputVariantIndex true [missing this arg. Needed to add extra step of generating index for genomicsDB import through IndexFeatureFile utility]
SPARK="-- --spark-runner LOCAL --spark-master local[$CPUS]"
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
        CMD="$CMD --intervals="$BED" --interval-padding 0"
    else
        printf "Invalid bed input for target regions.\n"
        display_usage
        exit 1
    fi  
fi


##Dbsnp file###
cat > $WORKDIR/haplotypeCaller.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=gatk_singleSample
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
srun="srun -N1 -n1 --slurmd-debug 3 shifter"
parallel="parallel --delay 0.2 -j \$SLURM_NTASKS  --env timestamp --joblog parallel.log --resume-failed"
echo "\$(timestamp) -> Job started!"
echo "\$(timestamp) -> Converting reference fasta file to 2bit format for HaplotypeCallerSpark"
\$srun --image=mcfonsecalab/htstools_plus:latest faToTwoBit $REF ${ref_2bit%.*}.2bit
echo "\$(timestamp) -> Done"
echo "\$(timestamp) -> Running HaplotypeCaller in normal mode"
if [ "$PARALLEL" == "true" ]; then
    cat $BAM_DATA | \$parallel '\$srun $CMD -I={} -O={/.}.vcf $SPARK && echo -e "{/.}\\t{/.}.vcf" >> aux_sampleName.map && shifter gatk IndexFeatureFile -F {/.}.vcf'
#    cat $BAM_DATA | \$parallel '\$srun $CMD -I={} -O={/.}.vcf && shifter gatk IndexFeatureFile -F {/.}.vcf'
else
    unique_samples=()
    for j in \$(find $BAM_DATA -exec cat {} \; );do
        printf "\$(timestamp): Processing \$(basename \$j) file.\n"
        i=\$(basename \$j)
        out=\$(echo "\$i" | cut -f1 -d "_")
        if [[ " \${unique_samples[@]} " =~ " \${out} " ]]; then
            printf "\$(timestamp): Warning, duplicate sample ID (\${out}) found after splitting by first '_'. Will use the filename instead without the extension."
            out="\${i%.*}"
        fi
        unique_samples+=(\${out})
        \$srun $CMD -I=\$j -O=\${out}.vcf -L chr1 $SPARK
#        \$srun $CMD -I=\$j -O=\${out}.vcf
        \$srun gatk IndexFeatureFile -F \${out}.vcf
        printf "\$(timestamp): Done!\n"
    done
fi
echo -e "\$(timestamp) -> All done!"
EOL
sbatch $WORKDIR/haplotypeCaller.sbatch
sleep 1 

cd $WORKDIR
mv haplotypeCaller.sbatch $(ls -td -- */ | head -n 1)
