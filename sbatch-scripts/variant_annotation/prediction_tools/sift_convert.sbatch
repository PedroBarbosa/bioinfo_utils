#!/bin/bash
#SBATCH --job-name=sift_convert
#SBATCH --time=24:00:00
#SBATCH --mem=150G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --workdir=/home/pedro.barbosa/scratch
#SBATCH --output=/home/pedro.barbosa/scratch/%j_sift_convert.log

while read line; do outchr=$(basename $line | cut -f1 -d "."); zcat $line | awk '($11 != "NA")' | grep -v "^#" | awk '{print $1-1"\t"$1"\t"$2$3"\t"$11"\t"$5"_"$6"_"$4}' | sed "s|^|chr${outchr}\t|" | gzip > chr${outchr}.bed.gz; done < /home/pedro.barbosa/resources/predictionTools/sift/GRCh37.74/listSiftFiles.txt
