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

srun shifter STAR --runThreadN $SLURM_CPUS_PER_TASK --runMode genomeGenerate --genomeDir $PWD --genomeFastaFiles /home/pedro.barbosa/mcfonseca/shared/genomes/hg38/GRCh38.primary.genome.fa --limitGenomeGenerateRAM 214748364800 --sjdbGTFfile /home/pedro.barbosa/mcfonseca/shared/genomes/hg38/gencode.v28.annotation.gtf --sjdbGTFtagExonParentTranscript transcript_id
