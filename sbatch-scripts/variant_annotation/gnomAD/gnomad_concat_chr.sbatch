#!/bin/bash
#SBATCH --job-name=gnomAD_concat
#SBATCH --time=72:00:00
#SBATCH --mem=240G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
##SBATCH --workdir=/home/pedro.barbosa/scratch
#SBATCH --output=%j_gnomAD_concat.log
#SBATCH --image=mcfonsecalab/variantutils:0.4

files2merge=$(readlink -f $1)
outfile="gnomAD_v2.1_justImportantFields.vcf.gz"
cd /home/pedro.barbosa/scratch
srun shifter bcftools concat -a -f $files2merge -Oz -o $outfile

mv $outfile $(dirname $files2merge)
rm /tmp/*
