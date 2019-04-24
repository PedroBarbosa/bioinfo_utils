#!/bin/bash

#SBATCH --job-name=bamInsight
#SBATCH --array=0-12%1
#SBATCH --time=72:00:00
#SBATCH --mem=240G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
##SBATCH --workdir=/home/pedro.barbosa/scratch/bamInsight
#SBATCH --image=rluis/baminsight:latest 
#SBATCH --output=%j_bamInsight.log

bamFiles=$(readlink -f $1) 
genome=$2
outdir=$3
readarray -t bams < $bamFiles

cd /home/pedro.barbosa/scratch/bamInsight
mkdir $SLURM_JOB_ID && cd $SLURM_JOB_ID

srun shifter baminsight stranded $genome pedro.barbosa@medicina.ulisboa.pt --names ${bams[$SLURM_ARRAY_TASK_ID]} --cpus $SLURM_CPUS_PER_TASK -FF 83 163 -FR 99 147 -FTPHOST immftp01.fm.ul.pt -FTPUser mcfonsecaftp -FTPPassword Lark699-evan -FTPPath /pedro.barbosa/$outdir/
# -FF 99 147 -FR 83 163 fr-secondstrand
#srun shifter baminsight stranded hg38 ruisergiolusi@gmail.com --names $file -FF 99 147 -FR 83 163 -FTPHOST immftp01.fm.ul.pt -FTPUser mcfonsecaftp -FTPPassword Lark699-evan -FTPPath $2  --cpus 40

#-add_bam
rm -rf /tmp/*
