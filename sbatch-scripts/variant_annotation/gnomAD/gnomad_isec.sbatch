#!/bin/bash
#SBATCH --job-name=gnomAD_isec
#SBATCH --time=72:00:00
#SBATCH --mem=240G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
##SBATCH --workdir=/home/pedro.barbosa/scratch
#SBATCH --output=%j_gnomAD_isec.log
#SBATCH --image=mcfonsecalab/variantutils:0.4

file1=$(readlink -f $1)
file2=$(readlink -f $2)
outdir=$(readlink -f $3)
#cd /home/pedro.barbosa/scratch
#zcat $file1 | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' | shifter --image=ummidock/ubuntu_base:latest bgzip > gnomAD_v2.1_justImportantFields_chrRenamed.vcf.gz
srun shifter bcftools isec -Oz -p $outdir -n=2 -w1 $file1 $file2

#mv $outfile $(dirname $files2merge)
rm /tmp/*
