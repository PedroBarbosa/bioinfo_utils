#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=7
#SBATCH --array=0-127%20
#SBATCH --mem=100G
#SBATCH --job-name="stringtie"
#SBATCH --output=%j_stringtie.out
#SBATCH --mail-user=pedro.barbosa@medicina.ulisboa.pt
#SBATCH --mail-type=ALL
#SBATCH --image=mcfonsecalab/hisat2pipeline:latest
##SBATCH --workdir=/home/pedro.barbosa/scratch/rna_seq/stringtie

listbams=$(readlink -f $1)
annotation=$(readlink -f $2)
out_directory=$(readlink -f $3)
mkdir /home/pedro.barbosa/scratch/rna_seq/stringtie/$SLURM_JOB_ID 
cd /home/pedro.barbosa/scratch/rna_seq/stringtie/$SLURM_JOB_ID

## --rf	Assumes a stranded library fr-firststrand.
## --fr	Assumes a stranded library fr-secondstrand.

readarray -t bams < <(cat $(readlink -f "$listbams"))
srun="srun -N1 -n1"
OUT_BASENAME=$(basename ${bams[$SLURM_ARRAY_TASK_ID]} | cut -f1,2,3,4,5 -d "_")
OUT_BASENAME=${OUT_BASENAME/Aligned.sortedByCoord.out.bam/}

srun shifter stringtie -p $SLURM_CPUS_PER_TASK -A ${OUT_BASENAME}_quantification.tsv -C ${OUT_BASENAME}_transcripts.tsv -b ${OUT_BASENAME}_ballgown -e -G $annotation -o ${OUT_BASENAME}.gtf ${bams[$SLURM_ARRAY_TASK_ID]}

#With parallel
#cat $listbams | parallel --tmpdir /home/pedro.barbosa/scratch/rna_seq/stringtie/ --delay 0.2 -j3 "srun shifter stringtie -p 15 -A {=s{.*/}{};s/\_[^_]+$//;=}_quantification.tsv -C {=s{.*/}{};s/\_[^_]+$//;=}_transcripts.tsv -b {=s{.*/}{};s/\_[^_]+$//;=}_ballgown -e -G $annotation -o {=s{.*/}{};s/\_[^_]+$//;=}.gtf {}"
mv * $out_directory
cd .. && rmdir $SLURM_JOB_ID
rm -rf `ls -la /tmp/ | grep 'pedro.barbosa' | awk ' { print $9 } '`
