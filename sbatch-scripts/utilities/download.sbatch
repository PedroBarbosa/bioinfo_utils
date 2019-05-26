#!/bin/bash
#SBATCH --job-name=download
#SBATCH --time=72:00:00
#SBATCH --mem=100Gb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --workdir=/home/pedro.barbosa/scratch/downloads
#SBATCH --output=/home/pedro.barbosa/scratch/%j_download.log

#wget ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad.genomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz
#wget ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad.genomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz.tbi
wget -c https://s3.eu-west-1.amazonaws.com/novogene-europe/HW/project/X201SC19042005-Z01-F001_20190522_ctepqH.tar
mv * $SLURM_SUBMIT_DIR 

