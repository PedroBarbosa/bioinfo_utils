#!/bin/bash
#SBATCH --job-name=generate_pop_gnomAD_genotypes
#SBATCH --array=0-2%3
#SBATCH --time=72:00:00
#SBATCH --mem=75G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=%j_generate_gnomAD_genotypes.log
#SBATCH --image=mcfonsecalab/variantutils:0.4

PRESENT_IN_GNOMAD=$(readlink -f $1)
ABSENT_IN_GNOMAD=$(readlink -f $2)
populations=$(readlink -f $3)
readarray -t pops < $populations
SCRIPT="/home/pedro.barbosa/git_repos/bioinfo_utils/python-scripts/vcf-tools/generate_gnomAD_genotypes.py"

cd /home/pedro.barbosa/scratch/vep
mkdir $SLURM_JOB_ID && cd $SLURM_JOB_ID

srun shifter python $SCRIPT --pop ${pops[$SLURM_ARRAY_TASK_ID]} $PRESENT_IN_GNOMAD $ABSENT_IN_GNOMAD gnomad_${pops[$SLURM_ARRAY_TASK_ID]}.vcf

srun shifter bcftools sort -Oz gnomad_${pops[$SLURM_ARRAY_TASK_ID]}.vcf > gnomad_${pops[$SLURM_ARRAY_TASK_ID]}.vcf.gz
#srun shifter --image=ummidock/ubuntu_base:latest bgzip --force gnomad_${pops[$SLURM_ARRAY_TASK_ID]}.vcf
srun shifter --image=ummidock/ubuntu_base:latest tabix -p vcf gnomad_${pops[$SLURM_ARRAY_TASK_ID]}.vcf.gz
mv gnomad_${pops[$SLURM_ARRAY_TASK_ID]}.vcf.gz* $(dirname $PRESENT_IN_GNOMAD)
rm -rf `ls -la /tmp/ | grep 'pedro.barbosa' | awk ' { print $9 } '`
