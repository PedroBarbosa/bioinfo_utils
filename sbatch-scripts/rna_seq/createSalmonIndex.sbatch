#!/bin/bash
#SBATCH --job-name=idx_salmon
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=%j_idx_salmon.log
#SBATCH --image=mcfonsecalab/rna_quantification:latest

fasta_file=$(readlink -f "$1")
idx_name=$2
srun shifter salmon index --gencode --keepDuplicates -i $idx_name -t $fasta_file

