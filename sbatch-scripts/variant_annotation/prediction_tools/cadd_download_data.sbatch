#!/bin/bash
#SBATCH --job-name=download_data
#SBATCH --time=24:00:00
#SBATCH --mem=220G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --workdir=/home/pedro.barbosa/scratch
#SBATCH --output=/home/pedro.barbosa/scratch/%j_cadd_download.log


#wget http://krishna.gs.washington.edu/download/CADD/v1.3/InDels_inclAnno.tsv.gz 
#wget http://krishna.gs.washington.edu/download/CADD/v1.3/InDels_inclAnno.tsv.gz.tbi
#wget http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz
#wget http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz.tbi
wget http://krishna.gs.washington.edu/download/CADD/v1.3/InDels.tsv.gz
wget http://krishna.gs.washington.edu/download/CADD/v1.3/InDels.tsv.gz.tbi
