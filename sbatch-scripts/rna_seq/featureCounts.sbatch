#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=100G
#SBATCH --job-name="featureCounts"
#SBATCH --output=%j_featureCounts.out
#SBATCH --mail-user=pedro.barbosa@medicina.ulisboa.pt
#SBATCH --image=mcfonsecalab/hisat2pipeline:latest
##SBATCH --workdir=/home/pedro.barbosa/scratch/rna_seq/stringtie

listbams=$(readlink -f $1)
annotation=$(readlink -f $2)
out_directory=$(readlink -f $3)
fasta=$(readlink -f "/home/pedro.barbosa/mcfonseca/shared/genomes/human/hg38/GRCh38.primary.genome.fa")
#fasta=$(readlink -f ~/mcfonseca/shared/genomes/human/hg38/all_scaffols_patches_haplotypes/GRCh38.p12.all_with_ERCC_spike-in.fa)
#fasta=$(readlink -f "/home/pedro.barbosa/mcfonseca/shared/genomes/mouse/GRCm38.p6/GRCm38.primary_assembly.genome.fa")
mkdir /home/pedro.barbosa/scratch/rna_seq/featureCounts/$SLURM_JOB_ID && cd /home/pedro.barbosa/scratch/rna_seq/featureCounts/$SLURM_JOB_ID

readarray -t bams < $listbams
files=""
for bam in "${bams[@]}"
do
    files+=" $bam"
done
#-t feature type: exon by default
#-J count reads supporting exon-exon juncion
#-G fasta of the reference sequences used in read mapping to improve read counting for junctions
#-s <int or string>  Perform strand-specific read counting. 0 (unstranded), 1 (stranded) and 2 (reversely stranded).
#--primary count primary alignments only identified using bit 0x100 in SAM/BAM FLAG field.
# --ignoreDup ignore duplicate reads identified using bit Ox400 in BAM/SAM FLAG field 
#--fracOverlap 0.2 --minOverlap 20
cmd="srun shifter featureCounts -T $SLURM_CPUS_PER_TASK --extraAttributes gene_name -p -s 2 --fracOverlap 0.1 --minOverlap 10 -J -G $fasta -Q 10 -a $annotation -o featureCounts.txt $files"
echo "-----COMMAND------"
echo "$cmd"
time $cmd
#time cat $listbams | parallel --tmpDir /home/pedro.barbosa/scratch/rna_seq/featureCounts/ --joblog parallel.log --delay 0.1 -j $SLURM_NTASKS_PER_NODE "srun shifter featureCounts -T $SLURM_CPUS_PER_TASK  --extraAttributes gene_name --fracOverlap 0.2 --minOverlap 20 -Q 1 --primary -a $annotation -o {/.}.txt {}"

mv * $out_directory
cd .. && rm -rf $SLURM_JOB_ID
rm -rf `ls -la /tmp/ | grep 'pedro.barbosa' | awk ' { print $9 } '`
