#!/bin/bash
#SBATCH --job-name=download
#SBATCH --time=72:00:00
#SBATCH --mem=100Gb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --workdir=/home/pedro.barbosa/scratch
#SBATCH --output=/home/pedro.barbosa/scratch/%j_download.log

wget -c https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.vcf.bgz
wget -c https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.vcf.bgz.tbi
mv gnomad.genomes.r2.1* mv ${SLURM_JOB_ID}_download.log /mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/gnomeAD/hg19

