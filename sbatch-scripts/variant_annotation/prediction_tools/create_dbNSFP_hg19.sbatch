#!/bin/bash
#SBATCH --job-name=concat
#SBATCH --time=24:00:00
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
##SBATCH --workdir=/homeo/pedro.barbosa/scratch
#SBATCH --output=%j_concat.sbatch.log

#header=$(readlink -f "$1")
export TMPDIR="/home/pedro.barbosa/scratch/"
file=$(readlink -f "$1")

srun unzip $file
head -n1 dbNSFP3.5a_variant.chr1 > h
cat dbNSFP3.5a_variant.chr* | grep -v "^#chr" | awk '$8 != "."' | sort -k8,8 -k9,9n - | cat h - | shifter --image=ummidock/ubuntu_base:latest bgzip -c > dbNSFP_hg19.gz
srun shifter --image=ummidock/ubuntu_base:latest tabix -s 8 -b 9 -e 9 dbNSFP_hg19.gz

#srun zcat $file | grep -v "^#" | awk '$8 != "."' | sort -k8,8 -k9,9n - | cat h - | shifter --image=ummidock/ubuntu_base:latest bgzip > dbNSFP_hg19.txt.gz
