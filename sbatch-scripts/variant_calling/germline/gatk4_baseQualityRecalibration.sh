#!/usr/bin/env bash
display_usage(){
 printf "Script to run GATK4 base quality score recalibration pipeline in lobo.\n
 Usage:.
    -1st argument must be list of bam files to recalibrate.
    -2nd argument must be the name of the final ouptut directory.
    -3th argument must be the fasta reference for which the reads were aligned. If '-' is set, the existing hg38 version in lobo will be used.
     Be aware that a fasta index file (.fai) must also be present in the fasta directory.
    -4th argument is optional. It refers to an additional resource of known sites to use in the calibration step. By default, hg38 1000Genomes and dbSNP SNPs are used."
 }


if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
    printf "Error. Please set the required parameters of the script\n"
    display_usage
    exit 1
fi

BAMS=$(readlink -f "$1")
OUTDIR=$(readlink -f "$2")
WORKDIR="/home/pedro.barbosa/scratch/gatk/BQSR"
if [ ! -d $WORKDIR ];then
    mkdir $WORKDIR
fi
if [ ! -d $OUTDIR ];then
    mkdir $OUTDIR
fi

##REFERENCE##
if [ "$3" = "-" ]; then
    REF="/home/pedro.barbosa/mcfonseca/shared/genomes/human/hg38/GRCh38.primary.genome.fa"
elif [ ! -f "$3" ]; then
    printf "Please provide a valid fasta file in th 4th argument.\n"
    display_usage
    exit 1
elif [ ! -f "${3}.fai" ]; then
    printf "Fasta index ${3}.fai not found in the reference directory. Please create one with samtools faidx.\n"
    display_usage
    exit
else
    REF="$3"
fi

known_sites="--known-sites /home/pedro.barbosa/resources/gatk-bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
--known-sites /home/pedro.barbosa/resources/gatk-bundle/hg38/dbsnp_146.hg38.vcf.gz "

if [[ ! -z "$4" ]]; then
    known_sites="$known_sites --known-sites $(readlink -f "$4")"
fi



CMD_RECALIBRATOR="gatk BaseRecalibrator -R $REF $known_sites"


cat > $WORKDIR/gatk4_baseRecalibration.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=BQSR
#SBATCH --time=72:00:00
#SBATCH --mem=240G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-0%1
#SBATCH --cpus-per-task=40
#SBATCH --image=broadinstitute/gatk:latest
##SBATCH --workdir=$WORKDIR
#SBATCH --output=%j_BQSR_gatk.log


timestamp() {
    date +"%Y-%m-%d  %T"
}
export -f timestamp
scratch_out=$WORKDIR/\$SLURM_JOB_ID
mkdir \$scratch_out
cd \$scratch_out
srun="srun -N1 -n1 --slurmd-debug 3 shifter"
echo "\$(timestamp) -> Recalibration job started!"

readarray -t bams < $BAMS

OUTBASENAME=\$(basename \${bams[\$SLURM_ARRAY_TASK_ID]})

\$srun $CMD_RECALIBRATOR -I \${bams[\$SLURM_ARRAY_TASK_ID]} -O \${OUTBASENAME/bam/bqsr.table} 

echo "\$(timestamp) -> Done! Apllying BQSR now."
CMD_APPLY_BQSR="gatk ApplyBQSR -bqsr  \${OUTBASENAME/bam/bqsr.table} -I \${bams[\$SLURM_ARRAY_TASK_ID]} -O \${OUTBASENAME/bam/bqsr.bam}"

\$srun \$CMD_APPLY_BQSR
echo -e "\$(timestamp) -> Finished job."
mv * $OUTDIR
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
echo -e "\$(timestamp) -> All done!"
EOL
sbatch $WORKDIR/gatk4_baseRecalibration.sbatch
sleep 1
cd $WORKDIR
mv gatk4_baseRecalibration.sbatch $(ls -td -- */ | head -n 1)
