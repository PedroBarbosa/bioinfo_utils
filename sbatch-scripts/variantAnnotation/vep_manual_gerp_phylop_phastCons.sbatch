#!/bin/bash
#SBATCH --job-name=vep_cache
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
##SBATCH --workdir=/home/pedro.barbosa/scratch/vep
#SBATCH --output=%j_vep_cache.log

cd /home/pedro.barbosa/scratch/vep
mkdir $SLURM_JOB_ID && cd $SLURM_JOB_ID
vcf_file=$(readlink -f $1) 

srun shifter -V=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.barbosa/resources/ensembl_vep:/media --image=ensemblorg/ensembl-vep:latest vep -i $vcf_file -o vep.out.vcf --cache --dir /media/cache --offline --vcf --sift b --polyphen b --force_overwrite -a GRCh37 --hgvs --fasta /media/cache/homo_sapiens/92_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz --canonical --variant_class --regulatory --gene_phenotype --hgvsg --numbers --check_existing --af --max_af --keep_csq --vcf_info_field ANN -custom /media/custom_data/gnomeAD/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz,gnomADg,vcf,exact,0,AF_NFE --custom /media/custom_data/gerp/All_hg19_RS.bw,GERP,bigwig --custom /media/custom_data/phastcons100/hg19.100way.phastCons.bw,phastCons,bigwig --custom /media/custom_data/phyloP100/hg19.100way.phyloP100way.bw,phyloP,bigwig  --plugin CADD,/media/custom_data/cadd/whole_genome_SNVs.tsv.gz,/media/custom_data/cadd/InDels.tsv.gz --plugin Gwava,region,/media/custom_data/gwava_scores/gwava_scores.bed.gz

#--max_af REport the highest allele frequency observed in any population from 1000 genomes, ESP or gnomAD. Not used by default
#--check_existing flag to identify known variants colocated with user input. VEP's known variant cache is derived from Ensembl's variation database and contains variants from dbSNP and other sources
#--plugin CADD fails 
#--fasta /mnt/nfs/lobo/IMM-NFS/genomes/hg19/Sequence/WholeGenomeFasta/genome.fa #fail index is not found
#--variant_class: output the Sequence Ontology variant class.
#--filter_common: this will exclude variants that have a co-located existing variant with global AF > 0.01 (1%) in the 1000 Genomes global population. 
#--freq_pop: 1KG_ALL or 1KG_EUR
#--freq_freq: Allele frequency
#--canonical: adds a flag to the output indicating if the reported transcript is the canonical transcript for the gene
#--numbers: report the exon or intron number that a variant falls in
#--plugins CADD, dbNSFP

mv ../${SLURM_JOB_ID}* /home/pedro.barbosa/hcm/variantAnnotation/gatk4_hg19/fromHardFiltering
