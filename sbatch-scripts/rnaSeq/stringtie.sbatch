#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=230G
#SBATCH --job-name="stringtie"
#SBATCH --output=%j_stringtie.out
#SBATCH --mail-user=pedro.barbosa@medicina.ulisboa.pt
#SBATCH --mail-type=ALL
#SBATCH --image=rluis/hisat2pipeline:1.0
##SBATCH --workdir=/home/pedro.barbosa/scratch/rna_seq/stringtie

listbams=$(readlink -f $1)
annotation=$(readlink -f $2)
out_directory=$(readlink -f $2)
mkdir /home/pedro.barbosa/scratch/rna_seq/stringtie/$SLURM_JOB_ID && cd /home/pedro.barbosa/scratch/rna_seq/stringtie/$SLURM_JOB_ID

## --rf	Assumes a stranded library fr-firststrand.
## --fr	Assumes a stranded library fr-secondstrand.

cat $listbams | parallel --tmpdir /home/pedro.barbosa/scratch/rna_seq/stringtie/ --delay 0.2 -j3 "srun shifter stringtie -p 15 -A {=s{.*/}{};s/\_[^_]+$//;=}_quantification.tsv -C {=s{.*/}{};s/\_[^_]+$//;=}_transcripts.tsv -b {=s{.*/}{};s/\_[^_]+$//;=}_ballgown -e -G $annotation -o {=s{.*/}{};s/\_[^_]+$//;=}.gtf {}"

mv * $out_directory
cd .. & rmdir $SLURM_JOB_ID
