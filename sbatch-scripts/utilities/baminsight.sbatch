#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=120G
#SBATCH --job-name="BamInsight"
#SBATCH --output=slurm_%j_BamInsight.out
#SBATCH --mail-user=rluis@medicina.ulisboa.pt
#SBATCH --mail-type=ALL
#SBATCH --image=docker:rluis/baminsight:latest


file=$(readlink -f $1)

cd /home/pedro.barbosa/scratch

mkdir baminsight_$SLURM_JOB_ID
cd baminsight_$SLURM_JOB_ID

#echo srun shifter baminsight starnded hg38 ruisergioluis@gmail.com --names $file -FTPHOST immftp01.fm.ul.pt -FTPUser mcfonsecaftp -FTPPassword Lark699-evan -FTPPath $2 --cpus 40

srun shifter baminsight stranded hg38  ruisergiolusi@gmail.com --names $file -FF 99 147 -FR 83 163 -FTPHOST immftp01.fm.ul.pt -FTPUser mcfonsecaftp -FTPPassword Lark699-evan -FTPPath $2  --cpus 40

echo srun shifter baminsight stranded hg38  ruisergiolusi@gmail.com --names $file -FF 32 -FR 16 -FTPHOST immftp01.fm.ul.pt -FTPUser mcfonsecaftp -FTPPassword Lark699-evan -FTPPath $2  --cpus 40 -long_label $3
	
#srun shifter baminsight stranded hg38 ruisergioluis@gmail.com --names $file -FTPHOST immftp01.fm.ul.pt -FTPUser mcfonsecaftp -FTPPassword Lark699-evan -FTPPath $2 --cpus 40

#srun shifter baminsight stranded hg38 ruisergiolusi@gmail.com --names $file -FF 16 -notFF -FR 16 -FTPHOST immftp01.fm.ul.pt -FTPUser mcfonsecaftp -FTPPassword Lark699-evan -FTPPath $2 --cpus 40 -long_label $3


cd ..
rm -r baminsight_$SLURM_JOB_ID
