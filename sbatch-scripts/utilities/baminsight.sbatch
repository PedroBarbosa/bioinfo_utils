#!/bin/bash
#SBATCH --job-name=bamInsight
#SBATCH --time=72:00:00
#SBATCH --mem=240G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
##SBATCH --workdir=/home/pedro.barbosa/scratch/bamInsight
#SBATCH --image=rluis/baminsight:latest 
#SBATCH --output=%j_bamInsight.log

bamFiles=$(readlink -f $1) 
genome=$2
outdir=$3

cd /home/pedro.barbosa/scratch/bamInsight
mkdir $SLURM_JOB_ID && cd $SLURM_JOB_ID

while read line; do
    echo "srun shifter baminsight stranded $genome pedro.barbosa@medicina.ulisboa.pt --names $line --cpus $SLURM_CPUS_PER_TASK -FF 83 163 -FR 99 147 -FTPHOST immftp01.fm.ul.pt -FTPUser mcfonsecaftp -FTPPassword Lark699-evan -FTPPath /pedro.barbosa/$outdir/"
done < $bamFiles
# -FF 99 147 -FR 83 163 fr-secondstrand

#-add_bam
rm -rf /tmp/*
