#!/bin/bash
#SBATCH --job-name=bamToBigWig
#SBATCH --array=0-2%4
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --workdir=/home/pedro.barbosa/scratch/bamInsight
#SBATCH --image=argrosso/htstools:0.2.1 
#SBATCH --output=%j_bam2bigwig.log

bamFiles=$(readlink -f $1) 
genomeChromSize=$(readlink -f $2)
outdir=$(readlink -f $3)
readarray -t bams < $bamFiles

cd /home/pedro.barbosa/scratch/bamUtilities
out=$(basename ${bams[$SLURM_ARRAY_TASK_ID]})
srun shifter genomeCoverageBed -ibam -split -bg -i ${bams[$SLURM_ARRAY_TASK_ID]} > ${out/.bam/.bg} #-g $genomeChromSize
srun shifter sortBed -i ${out/.bam/.bg} > ${out/.bam/.sort.bg}
srun shifter bedGraphToBigWig ${out/.bam/.sort.bg} $genomeChromSize ${out/.bam/.bw}

mv ${out/.bam/.bw} $outdir
rm ${out/.bam/.bg} ${out/.bam/.sort.bg}
