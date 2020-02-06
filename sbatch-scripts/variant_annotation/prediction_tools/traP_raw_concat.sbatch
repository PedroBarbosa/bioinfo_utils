#!/bin/bash
#SBATCH --job-name=trap_convert
#SBATCH --time=72:00:00
#SBATCH --mem=150G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --output=/home/pedro.barbosa/scratch/%j_trap_convert.log
#SBATCH --image=mcfonsecalab/variantutils:0.4

chr_indiv_scores=$(readlink -f $1)
output_file=$(readlink -f $2)
cd /home/pedro.barbosa/scratch

while read line; do
    chr_name=$(echo $(basename $line) | cut -f1 -d ".")
    zcat $line | awk -v chr="$chr_name" '{ print chr,$0 }' OFS='\t'  >> $output_file
done < $chr_indiv_scores
gzip $output_file
mv ${SLURM_JOB_ID}* $(dirname $output_file)
