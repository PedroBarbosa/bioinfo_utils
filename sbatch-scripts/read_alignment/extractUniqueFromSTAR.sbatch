#!/bin/bash
#SBATCH --job-name=extractUnique
#SBATCH --array=0-8%5
#SBATCH --time=72:00:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
##SBATCH --workdir=/home/pedro.barbosa/scratch
#SBATCH --output=%j_extractUnique.log
#SBATCH --image=mcfonsecalab/htstools_plus:latest

bamFiles=$(readlink -f "$1")
outdir=$(readlink -f "$2")
readarray -t bams < $bamFiles

cd /home/pedro.barbosa/scratch

outbasename=$(basename ${bams[$SLURM_ARRAY_TASK_ID]})
srun shifter --image=mcfonsecalab/htstools_plus:latest samtools view -h --threads $SLURM_CPUS_PER_TASK ${bams[$SLURM_ARRAY_TASK_ID]} | grep -e "NH:i:1" -e "^@" | shifter samtools view -hb - >  ${outbasename/.bam/_uniq.bam}
mv ${outbasename/.bam/_uniq.bam} $outdir
find /tmp/ -user pedro.barbosa -exec rm -rf {} \;

