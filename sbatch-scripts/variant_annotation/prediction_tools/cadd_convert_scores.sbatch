#!/bin/bash
#SBATCH --job-name=cadd_convert
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --workdir=/home/pedro.barbosa/scratch
#SBATCH --output=/home/pedro.barbosa/scratch/%j_cadd_convert.log

filein=$(readlink -f /home/pedro.barbosa/resources/predictionTools/cadd13/whole_genome_SNVs.tsv.gz)
zcat $filein | grep -v "#" | awk 'OFS="\t" {print "chr"$1,$2-1,$2,$3 $4,$6}' | gzip > cadd13_whole_genome_SNVs.bed.gz
mv cadd13_whole_genome_SNVs.bed.gz /home/pedro.barbosa/resources/predictionTools/cadd13
rm ${SLURM_JOB_ID}_cadd_convert.log
