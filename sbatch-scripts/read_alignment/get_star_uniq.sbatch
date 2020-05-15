#!/bin/bash
#SBATCH --job-name=get_star_uniq
#SBATCH --time=24:00:00
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
##SBATCH --workdir=/home/pedro.barbosa/scratch
#SBATCH --output=uniq/%j_get_star_uniq.log

if [[ ! -d "$PWD/uniq" ]]; then
    mkdir uniq
fi
bam=$(readlink -f "$1") 

cd /home/pedro.barbosa/scratch
samtools view -hb -q255 --threads $SLURM_CPUS_PER_TASK $bam -o ${bam/.bam/_uniq.bam}
mv ${bam/.bam/_uniq.bam} $SLURM_SUBMIT_DIR/uniq
