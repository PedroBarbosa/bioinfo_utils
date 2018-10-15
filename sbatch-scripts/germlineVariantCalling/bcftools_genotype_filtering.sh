#!/usr/bin/env bash
display_usage(){
 printf "Script to run bcftools in parallel to filter individual genotypes.\n
 Usage:.
    -1st argument must be the input VCF file to process.
    -2nd argument must be the name of the final ouptut directory.
    -3rd argument musth be the basename of the output file.
    -4th argument is optional. It refers to a file including the samples to process. Default: Process all samples.\n"
 }

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] ; then
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
    printf "Unknown format of the VCF file. Must end with .vcf or .gz\n"
    display_usage
    exit 1
fi

####SCRATCH WORKDIR####
OUTDIR=$(readlink -f "$2")
OUTBASENAME="$3"

WORKDIR="/home/pedro.barbosa/scratch/gatk/bcftools_hardFilt"
if [ ! -d $WORKDIR ];then
    mkdir $WORKDIR
fi
if [ ! -d $OUTDIR ];then
    mkdir $OUTDIR
fi

#BCFTOOLS: bcftools filter -s LowGQ,LowAD -m + -S . -e "FMT/GQ < 50 " -m + -e "FMT/DP[1] <= 15" ~/Desktop/joint.vcf.gz

samples_subset=$(readlink -f "$4")
if [ -z "$4" ]; then
    shifter --image=mcfonsecalab/variantutils:0.4 bcftools query -l $VCF > $WORKDIR/listSamples.txt
elif [ -f "$samples_subset" ]; then
    cat $samples_subset > $WORKDIR/listSamples.txt
else
    printf "Please set a valid file for the 4th argument"
    display_usage
    exit 1
fi

JOBS=$(( $(cat $WORKDIR/listSamples.txt | wc -l) - 1 ))


cat > $WORKDIR/bcftools_genotypes_hardFiltering.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=filt_hard_bcftools
#SBATCH --array=0-$JOBS%10
#SBATCH --time=24:00:00
#SBATCH --mem=15G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --image=mcfonsecalab/variantutils:0.4
#SBATCH --workdir=$WORKDIR
#SBATCH --output=$WORKDIR/%j_hardFiltbcftools.log


scratch_out=$WORKDIR/\$SLURM_JOB_ID
mkdir \$scratch_out
cd \$scratch_out
srun="srun -N1 -n1"

readarray -t samples < $WORKDIR/listSamples.txt

\$srun shifter bcftools view -Oz -s \${samples[\$SLURM_ARRAY_TASK_ID]} -o \${samples[\$SLURM_ARRAY_TASK_ID]}.vcf.gz $VCF
\$srun shifter bcftools view --exclude-uncalled -i 'FMT/DP >= 20 & FMT/GQ > 30 & (FMT/GT="0/0" & FMT/AD[0:0] >= 20 || MIN(FMT/AD) > 7)' -Oz -o \${samples[\$SLURM_ARRAY_TASK_ID]}_filt.vcf.gz \${samples[\$SLURM_ARRAY_TASK_ID]}.vcf.gz
 \$srun shifter bcftools index \${samples[\$SLURM_ARRAY_TASK_ID]}_filt.vcf.gz

mv *_filt.vcf.gz* ../

echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
EOL
job_submission=$(sbatch $WORKDIR/bcftools_genotypes_hardFiltering.sbatch)
id=${job_submission##* }

cd $WORKDIR
echo -e "#!/bin/bash
#SBATCH --job-name=mergeVCFs
#SBATCH --time=72:00:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=mcfonsecalab/variantutils:0.4
#SBATCH --workdir=$WORKDIR
#SBATCH --output=$WORKDIR/merge_bcftools_%j.log\n" > $WORKDIR/mergeVCFs.sbatch

echo -e "srun shifter bcftools merge -Oz -i DP:min *_filt.vcf.gz | shifter bcftools norm -Oz -m -both -o ${OUTBASENAME}_merged.vcf.gz -\n" >> $WORKDIR/mergeVCFs.sbatch
echo -e "srun shifter bcftools index ${OUTBASENAME}_merged.vcf.gz\n" >> $WORKDIR/mergeVCFs.sbatch
echo -e "srun shifter bcftools view --min-ac 1 -Oz -o ${OUTBASENAME}_merged_filt.vcf.gz ${OUTBASENAME}_merged.vcf.gz\n" >> $WORKDIR/mergeVCFs.sbatch
echo -e "mv ${OUTBASENAME}_merged* merge*log $OUTDIR" >> $WORKDIR/mergeVCFs.sbatch
echo "find . ! -name 'mergeVCFs.sbatch' -exec rm -rf {} \;" >> $WORKDIR/mergeVCFs.sbatch
sbatch --depend=afterok:$id $WORKDIR/mergeVCFs.sbatch

