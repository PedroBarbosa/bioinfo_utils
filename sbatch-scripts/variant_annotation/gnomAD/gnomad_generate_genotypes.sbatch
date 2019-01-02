#!/bin/bash
#SBATCH --job-name=generate_pop_gnomAD_genotypes
#SBATCH --array=0-2%3
#SBATCH --time=72:00:00
#SBATCH --mem=150G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=%j_generate_gnomAD_genotypes.log
#SBATCH --image=mcfonsecalab/variantutils:0.4

PRESENT_IN_GNOMAD=$(readlink -f $1)
ABSENT_IN_GNOMAD=$(readlink -f $2)
populations=$(readlink -f $3)
original_vcf=$(readlink -f $4)
readarray -t pops < $populations
SCRIPT="/home/pedro.barbosa/git_repos/bioinfo_utils/python-scripts/vcf-tools/generate_gnomAD_genotypes.py"

cd /home/pedro.barbosa/scratch/vep
mkdir $SLURM_JOB_ID && cd $SLURM_JOB_ID

srun shifter python $SCRIPT --pop ${pops[$SLURM_ARRAY_TASK_ID]} $PRESENT_IN_GNOMAD $ABSENT_IN_GNOMAD gnomad_${pops[$SLURM_ARRAY_TASK_ID]}.vcf
srun shifter --image=ummidock/ubuntu_base:latest bgzip gnomad_${pops[$SLURM_ARRAY_TASK_ID]}.vcf && srun shifter --image=ummidock/ubuntu_base:latest tabix -p vcf gnomad_${pops[$SLURM_ARRAY_TASK_ID]}.vcf.gz
srun shifter bcftools merge -0 --merge all $original_vcf gnomad_${pops[$SLURM_ARRAY_TASK_ID]}.vcf.gz | shifter bcftools sort -Oz -o gnomad_merged_${pops[$SLURM_ARRAY_TASK_ID]}.vcf.gz
srun shifter --image=ummidock/ubuntu_base:latest tabix -p vcf gnomad_merged_${pops[$SLURM_ARRAY_TASK_ID]}.vcf.gz
mv gnomad_merged_${pops[$SLURM_ARRAY_TASK_ID]}.vcf.gz* $(dirname $PRESENT_IN_GNOMAD)
cd .. && rm -rf $SLURM_JOB_ID
rm -rf `ls -la /tmp/ | grep 'pedro.barbosa' | awk ' { print $9 } '`
