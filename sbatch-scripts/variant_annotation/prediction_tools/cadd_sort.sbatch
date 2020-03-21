#!/bin/bash
#SBATCH --job-name=cadd_sort
#SBATCH --time=24:00:00
#SBATCH --mem=240G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --workdir=/home/pedro.barbosa/scratch
#SBATCH --output=/home/pedro.barbosa/scratch/%j_cadd_sort.log

#filein=$(readlink -f /home/pedro.barbosa/resources/predictionTools/cadd13/cadd13_whole_genome_SNVs.bed.gz)
#zcat $filein | sort -T /home/pedro.barbosa/scratch -k1,1 -k2,2n  > cadd13_wgs_sorted_SNVs.bed.gz
#mv cadd13_wgs_sorted_SNVs.bed.gz /home/pedro.barbosa/resources/predictionTools/cadd13
gzip < /home/pedro.barbosa/resources/predictionTools/cadd13/cadd13_wgs_sorted_SNVs.bed > cadd13_wgs_sorted_SNVs.bed.gz 
mv cadd13_wgs_sorted_SNVs.bed.gz /home/pedro.barbosa/resources/predictionTools/cadd13/
mv ${SLURM_JOB_ID}_cadd_sort.log /home/pedro.barbosa/resources/predictionTools/cadd13/
