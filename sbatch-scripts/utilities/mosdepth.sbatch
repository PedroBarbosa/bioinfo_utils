#!/bin/bash
#SBATCH --job-name=mosdepth
#SBATCH --array=0-4%5
#SBATCH --time=24:00:00
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
##SBATCH --workdir=/home/pedro.barbosa/scratch/covAnalysis
#SBATCH --output=%j_mosdepth.log
#SBATCH --image=mcfonsecalab/coverage_tools:latest 

bams=$(readlink -f $1)
outdir=$(readlink -f $2)
regions_file=$(readlink -f $3)
readarray -t files < $bams

#"${filename%.*}"
echo $(echo ${files[$SLURM_ARRAY_TASK_ID]##*/} | cut -f1 -d ".")
echo ${files[$SLURM_ARRAY_TASK_ID]} 
cd /home/pedro.barbosa/scratch/covAnalysis && mkdir $SLURM_JOB_ID && cd $SLURM_JOB_ID

srun shifter mosdepth --mapq 1 --thresholds 20,50,75,90,95 --by $regions_file -n $(echo ${files[$SLURM_ARRAY_TASK_ID]##*/} | cut -f1 -d ".") ${files[$SLURM_ARRAY_TASK_ID]}


mv * $outdir && cd ../ && rm -rf $SLURM_JOB_ID

