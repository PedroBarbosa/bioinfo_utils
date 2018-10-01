#!/usr/bin/env bash
display_usage(){
 printf "Script to run GATK4 Hard Filtering procedures.\n
 Usage:.
    -1st argument must be the input VCF file to recalibrate.
    -2nd argument must be the name of the final ouptut directory.
    -3rd argument musth be the name of the output file.
    -4th argument must the the mode to apply the filtering: (SNP,INDEL,BOTH).
    -5th argument is optional. Refers to apply genotype based filters on the VCF. If true, It needs a list of samples in the 6th argument (e.g bcftools query -l vcfile).\
    Default: false. Values: [true|false]
    -6th argument is optional. May only be set if 5th argument is set to true. Refers to the list of samples to apply genotype-based filters.\n"
 }

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ]; then
    printf "Error. Please set the required parameters of the script\n"
    display_usage
    exit 1
fi

if [[ ${1: -4} == ".vcf" || ${1: -3} == ".gz" || ${1: -4} == ".bgz" ]] ; then
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

if [[ -z "$5" || "$5" = "false" ]]; then
    FORMAT_FILTER=""
elif [[ "$5" = "true" ]]; then
    if [[ -z "$6" || ! -f "$6" ]] ; then
        printf "When 5th argument is set, you need to provide a valid file in the 6th argument with the sample names to process (one per line).\n"
        display_usage
        exit 1
    else
        FORMAT_FILTER="--set-filtered-genotype-to-no-call true --genotype-filter-expression "GQ <= 20" --genotype-filter-name "LowGQ""
        while read line; do
            FORMAT_FILTER = "$FORMAT_FILTER --genotype-filter-expression \"vc.getGenotype(\"$line\").getAD().1 < 10 --genotype-filter-name \"${line}_LowAD\" \
            --genotype-filter-expression \"vc.getGenotype(\"$line\").getDP() < 30 --genotype-filter-name \"${line}_LowDP\" "
        done < "$6"
    fi
#BCFTOOLS: bcftools filter -s LowGQ,LowAD -m + -S . -e "FMT/GQ < 50 " -m + -e "FMT/DP[1] <= 15" ~/Desktop/joint.vcf.gz
#docker run -v /home/pbarbosa/Desktop/:/media broadinstitute/gatk:latest gatk SelectVariants --exclude-non-variants -V /media/test -O /media/test2.vcf
#docker run -v /home/pbarbosa/Desktop/:/media broadinstitute/gatk:latest gatk SelectVariants --set-filtered-gt-to-nocall true \
# --select "vc.getGenotype("dAgo1_3_6").getDP() <= 30" -V /media/test -O /media/test2.vcf

####SCRATCH WORKDIR####
OUTDIR=$(readlink -f "$2")
OUTFILE="$3"
WORKDIR="/home/pedro.barbosa/scratch/gatk/hardFilt"
if [ ! -d $WORKDIR ];then
    mkdir $WORKDIR
fi
if [ ! -d $OUTDIR ];then
    mkdir $OUTDIR
fi

##MODE
MODE_ARRAY=("SNP" "INDEL" "BOTH")
if [[ !  " ${MODE_ARRAY[@]} " =~ " ${4} " ]]; then
    printf "Please set a valid value for the 4th argument.\n"
    display_usage
    exit 1
else
    MODE="$4"
fi

cat > $WORKDIR/gatk4_hardFiltering.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=filt_hard_gatk
#SBATCH --time=72:00:00
#SBATCH --mem=240G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --image=broadinstitute/gatk:latest
#SBATCH --workdir=$WORKDIR
#SBATCH --output=$WORKDIR/%j_hardFiltgatk.log


scratch_out=$WORKDIR/\$SLURM_JOB_ID
mkdir \$scratch_out
cd \$scratch_out
srun="srun -N1 -n1"
if [[ "$MODE" = "SNP" || "$MODE" = "BOTH" ]]; then

    \$srun shifter gatk VariantFiltration --filter-name "LowQD" --filter-expression "QD < 2.0"  --filter-name LowMQ --filter-expression "MQ < 40.0" --filter-name HighFS --filter-expression "FS > 60.0" --filter-name LowMQRankSum --filter-expression "MQRankSum < -12.5" --filter-name HighSOR --filter-expression "SOR > 3.0" --filter-name LowReadPosRankSum --filter-expression "ReadPosRankSum < -8.0" -O $OUTFILE -V $VCF
#    mv $OUTFILE $OUTDIR $FORMAT_FILTER
elif [[ "$MODE" = "INDEL" ]]; then
    \$srun shifter gatk VariantFiltration --filter-name LowQD --filter-expression "QD < 2.0" --filter-name LowMQ --filter-expression "MQ < 40.0" --filter-name HighFS --filter-expression "FS > 200.0" --filter-name LowMQRankSum --filter-expression "MQRankSum < -12.5" --filter-name HighSOR --filter-expression "SOR > 10.0" --filter-name LowReadPosRankSum --filter-expression "ReadPosRankSum < -20.0" -O $OUTFILE -V $VCF
#    mv $OUTFILE $OUTDIR $FORMAT_FILTER
fi
mv * ../\${SLURM_JOB_ID}_hardFiltgatk.log $OUTDIR

echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
EOL
sbatch $WORKDIR/gatk4_hardFiltering.sbatch
sleep 1
cd $WORKDIR
mv gatk4_hardFiltering.sbatch $(ls -td -- */ | head -n 1)

