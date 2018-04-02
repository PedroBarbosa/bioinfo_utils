#!/usr/bin/env bash
display_usage(){
 printf "Script to run GATK4 Variant Recalibration pipeline in lobo.\n
 Usage:.
    -1st argument must be the input VCF file to recalibrate.
    -2nd argument must be the name of the final ouptut directory.
    -3rd argument must refer to the type of analysis in hand. Values:[WGS|targeted]. Useful to apply different values in some arguments.
    -4th argument must be the fasta reference for which the reads were aligned. If '-' is set, the existing hg38 version in lobo will be used.
     Be aware that a fasta index file (.fai) must also be present in the fasta directory.
    -5th argument must be a file listing the resources to employ for in the recalibration. Must be a 2 column tab separated file where first
    column must refer to the input file and second to the resource name. Available resources names: [hapmap,omni,1000G,dbsnp,mills]. At least
     one resource is mandatory.
    -6th argument must the the mode to apply recalibration: (SNP,INDEL,BOTH).\n"
 }

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ]; then
    printf "Error. Please set the required parameters of the script\n"
    display_usage
    exit 1
fi

if [[ ${1: -4} == ".vcf" || ${1: -3} == ".gz" ]]; then
    if [ -f "$1" ]; then
        VCF=$(readlink -f "$1")
    else
        printf "VCF file not found!! ($1)"
        display_usage
        exit 1
    fi
else
    printf "Unknown format of the VCF file. Must end with .vcf or .gz"
    display_usage
    exit 1
fi

####SCRATCH WORKDIR####
OUTDIR=$(readlink -f "$2")
WORKDIR="/home/pedro.barbosa/scratch/gatk/VSQR"
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

##EXPERIMENT
analysis=("WGS" "targeted")
if [[ !  " ${analysis[@]} " =~ " ${3} " ]]; then
    printf "Please set a valid value for the 3th argument.\n"
    display_usage
    exit 1
else
    experiment="$3"
fi

CMD_RECALIBRATOR="gatk VariantRecalibrator -R $REF -V $VCF"
##MODE
MODE_ARRAY=("SNP" "INDEL" "BOTH")
if [[ !  " ${MODE_ARRAY[@]} " =~ " ${6} " ]]; then
    printf "Please set a valid value for the 6th argument.\n"
    display_usage
    exit 1
else
    MODE="$6"
fi
##RESOURCES
if [ -f "$5" ]; then
    while read line; do
        resource_file=$(echo "$line" | cut -f1)
        if [ ! -f $resource_file ]; then
            printf "$resource_file provided in the $5 argument does not exist."
            display_usage
            exit 1
        fi
        db=$(echo "$line" | cut -f2)
        if [ $MODE == "SNP" ] || [ $MODE == "BOTH" ]; then
            case "$db" in
                "hapmap")
                    CMD_RECALIBRATOR="$CMD_RECALIBRATOR -resource hapmap,known=false,training=true,truth=true,prior=15.0:$resource_file"
                    ;;
                "omni")
                    CMD_RECALIBRATOR="$CMD_RECALIBRATOR -resource omni,known=false,training=true,truth=true,prior=12.0:$resource_file"
                    ;;
                "dbsnp")
                    CMD_RECALIBRATOR="$CMD_RECALIBRATOR -resource 1000G,known=false,training=true,truth=false,prior=10.0:$resource_file"
                    ;;
                "1000G")
                    CMD_RECALIBRATOR="$CMD_RECALIBRATOR -resource dbsnp,known=true,training=false,truth=false,prior=2.0:$resource_file"
                    ;;
                *)
                    printf "Invalid database provided in the list of resources"
                    display_usage
                    exit 1
                    ;;
            esac
        elif [ $MODE = "INDEL" ]; then
            case "$db" in
                "mills")
                    CMD_RECALIBRATOR="$CMD_RECALIBRATOR -resource mills,known=false,training=true,truth=true,prior=12.0:$resource_file"
                    ;;
                "dbsnp")
                    CMD_RECALIBRATOR="$CMD_RECALIBRATOR -resource dbsnp,known=false,training=true,truth=false,prior=10.0:$resource_file"
                    ;;
                *)
                    printf "Invalid database provided in the list of resources"
                    display_usage
                    exit 1
                    ;;
            esac
        fi
    done < "$5"
else
    printf "$5 file does not exist."
    display_usage
    exit 1
fi

##ANNOTATIONS
coverage="-an DP" #Total depth of coverage
qualbydepth="-an QD"  #Variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples.
fisherstrand="-an FS" #This is the Phred-scaled probability that there is strand bias at the site
strandOddRatio="-an SOR" #Another way to test for strand bias. FS tends to penalize variants that occur at the ends of exons.
rmsMappingQuality="-an MQ" #Square root of the average of the squares of the mapping qualities at the site. It is meant to
# include the standard deviation of the mapping qualities.
mappingQualRankSumTest="-an MQRankSum" #It compares the mapping qualities of the reads supporting the reference allele and the alternate allele.
readPosRankSumTest="-an ReadPosRankSum" # It compares whether the positions of the reference and alternate alleles are different within the reads. Seeing an allele only near the ends
# of reads is indicative of error, because that is where sequencers tend to make the most errors.

CMD_RECALIBRATOR="$CMD_RECALIBRATOR $qualbydepth $fisherstrand $strandOddRatio $rmsMappingQuality $mappingQualRankSumTest $readPosRankSumTest"
if [ "$experiment" = "WGS" ]; then
    CMD_RECALIBRATOR="$CMD_RECALIBRATOR $coverage"
else
    CMD_RECALIBRATOR="$CMD_RECALIBRATOR --max-gaussians 4" 
fi


##TRANCHES
#Tranches are essentially slices of variants, ranked by VQSLOD, bounded by the threshold values specified
# in this step. The threshold values themselves refer to the sensitivity we can obtain when we apply them
# to the call sets that the program uses to train the model.
SENSITIVITY_THRESHOLDS="-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0"
CMD_RECALIBRATOR="$CMD_RECALIBRATOR $SENSITIVITY_THRESHOLDS"

cat > $WORKDIR/gatk4_recalibration.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=VQSR
#SBATCH --time=72:00:00
#SBATCH --mem=240G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --image=broadinstitute/gatk:latest
#SBATCH --workdir=$WORKDIR
#SBATCH --output=$WORKDIR/%j_VQSR_gatk.log


timestamp() {
    date +"%Y-%m-%d  %T"
}
export -f timestamp
scratch_out=$WORKDIR/\$SLURM_JOB_ID
mkdir \$scratch_out
cd \$scratch_out
srun="srun -N1 -n1 --slurmd-debug 3 shifter"
echo "\$(timestamp) -> Recalibration job started!"
\$srun $CMD_RECALIBRATOR -O recalibration.out -tranches-file tranches.out --rscript-file output.plots.R
echo "\$(timestamp) -> Done! Apllying VQSR now."
CMD_APPLY_VQSR="gatk ApplyVQSR -R $REF -V $VCF --recal-file recalibration.out --tranches-file tranches.out -O out_vqsr.vcf.gz"

if [ "$MODE" = "SNP" ]; then
    CMD_APPLY_VQSR="\$CMD_APPLY_VQSR --truth-sensitivity-filter-level 99.5 -mode SNP"
elif [ "$MODE" = "INDEL" ]; then
    CMD_APPLY_VQSR="\$CMD_APPLY_VQSR --truth-sensitivity-filter-level 90.0 -mode INDEL"
else
    CMD_APPLY_VQSR="\$CMD_APPLY_VQSR --truth-sensitivity-filter-level 95.0 -mode BOTH"
fi
\$srun \$CMD_APPLY_VQSR
echo -e "\$(timestamp) -> Finished job."
mv * ../\${SLURM_JOB_ID}_VQSR_gatk.log $OUTDIR
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
echo -e "\$(timestamp) -> All done!"
EOL
sbatch $WORKDIR/gatk4_recalibration.sbatch
sleep 1
cd $WORKDIR
mv gatk4_recalibration.sbatch $(ls -td -- */ | head -n 1)
