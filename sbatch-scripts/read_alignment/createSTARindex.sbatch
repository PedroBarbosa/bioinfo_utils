#!/bin/bash
#SBATCH --job-name=createStartIdx
#SBATCH --time=24:00:00
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
##SBATCH --workdir=/home/pedro.barbosa/scratch
#SBATCH --output=%j_createStartIdx.log
#SBATCH --image=argrosso/star:2.6.0b

fasta=$(readlink -f "$1")
annotation=$(readlink -f "$2")
srun shifter STAR --runThreadN $SLURM_CPUS_PER_TASK --runMode genomeGenerate --genomeDir $PWD --genomeFastaFiles $fasta --limitGenomeGenerateRAM 214748364800 --sjdbGTFfile $annotation --sjdbGTFtagExonParentTranscript transcript_id
