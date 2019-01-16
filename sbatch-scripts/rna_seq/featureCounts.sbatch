#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=230G
#SBATCH --job-name="featureCounts"
#SBATCH --output=%j_featureCounts.out
#SBATCH --mail-user=pedro.barbosa@medicina.ulisboa.pt
#SBATCH --image=mcfonsecalab/hisat2pipeline:latest
##SBATCH --workdir=/home/pedro.barbosa/scratch/rna_seq/stringtie

listbams=$(readlink -f $1)
annotation=$(readlink -f $2)
out_directory=$(readlink -f $3)
mkdir /home/pedro.barbosa/scratch/rna_seq/featureCounts/$SLURM_JOB_ID && cd /home/pedro.barbosa/scratch/rna_seq/featureCounts/$SLURM_JOB_ID

readarray -t bams < $listbams
files=""
for bam in "${bams[@]}"
do
    files+=" $bam"
done
#-t feature type: exon by default
time srun shifter featureCounts -T $SLURM_CPUS_PER_TASK --extraAttributes gene_name --fracOverlap 0.2 --minOverlap 20 -Q 1 --primary -a $annotation -o featureCounts.txt $files
#time cat $listbams | parallel --tmpDir /home/pedro.barbosa/scratch/rna_seq/featureCounts/ --joblog parallel.log --delay 0.1 -j $SLURM_NTASKS_PER_NODE "srun shifter featureCounts -T $SLURM_CPUS_PER_TASK  --extraAttributes gene_name --fracOverlap 0.2 --minOverlap 20 -Q 1 --primary -a $annotation -o {/.}.txt {}"

mv * $out_directory
cd .. && rm -rf $SLURM_JOB_ID
rm -rf `ls -la /tmp/ | grep 'pedro.barbosa' | awk ' { print $9 } '`
