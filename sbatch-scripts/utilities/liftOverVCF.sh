#!/usr/bin/env bash
display_usage(){
 printf "Script to run liftOver variant tools from GATK4 toolkit.\n
 Usage:.
    -1st argument must be the list of input VCF files to liftover.
    -2nd argument must be the name of the final ouptut directory.
    -3rd argument must be the fasta reference sequence.
    -4th argument must be the chain file.\n" 
}

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ]; then
    printf "Error. Please set the required parameters of the script\n"
    display_usage
    exit 1
fi


####SCRATCH WORKDIR####
VCFs=$(readlink -f "$1")
OUTDIR=$(readlink -f "$2")
REF=$(readlink -f "$3")
CHAIN=$(readlink -f "$4")

WORKDIR="/home/pedro.barbosa/scratch/gatk/"
if [ ! -d $WORKDIR ];then
    mkdir $WORKDIR
fi
if [ ! -d $OUTDIR ];then
    mkdir $OUTDIR
fi

CMD="gatk LiftoverVcf --CHAIN $CHAIN --REFERENCE_SEQUENCE $REF"
cat > $WORKDIR/gatk4_liftOverVCFs.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=vcfLiftOver
#SBATCH --time=72:00:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=8
#SBATCH --image=broadinstitute/gatk:latest
#SBATCH --workdir=$WORKDIR
#SBATCH --output=$WORKDIR/%j_vcfLiftOver.log


scratch_out=$WORKDIR/\$SLURM_JOB_ID
mkdir \$scratch_out
cd \$scratch_out
srun="srun -N1 -n1"

cat $VCFs | parallel --tmpdir $WORKDIR -j \$SLURM_NTASKS --joblog parallel.log --resume-failed "\$srun shifter $CMD --REJECT {/.}.rejected.vcf --INPUT {} --OUTPUT {/.}_liftOver.vcf"

mv * ../\${SLURM_JOB_ID}_liftOver.log $OUTDIR
cd ../ && rm -rf \$scratch_out
EOL
sbatch $WORKDIR/gatk4_liftOverVCFs.sbatch
sleep 5
cd $WORKDIR
mv gatk4_liftOverVCFs.sbatch $(ls -td -- */ | head -n 1)

