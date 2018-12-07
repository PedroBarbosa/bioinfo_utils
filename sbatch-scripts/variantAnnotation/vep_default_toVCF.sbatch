#!/bin/bash
#SBATCH --job-name=defaultFormat2VCF
#SBATCH --time=72:00:00
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
##SBATCH --workdir=/home/pedro.barbosa/scratch/vep
#SBATCH --output=%j_vep_cache.log

variants=$(readlink -f $1) 
output_vcf=$2
outdir=$(readlink -f $3)

cd /home/pedro.barbosa/scratch/vep && mkdir $SLURM_JOB_ID && cd $SLURM_JOB_ID

srun shifter -V=/mnt/nfs/lobo/IMM-NFS/ensembl_vep:/media --image=ensemblorg/ensembl-vep:latest vep -i $variants -o $output_vcf --compress_output bgzip --stats_text --hgvs --hgvsg --pubmed --cache --dir /media/cache --offline --vcf --numbers --regulatory  --sift b --polyphen b --force_overwrite -a GRCh37 --custom /media/custom_data/gnomeAD/hg19/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz,gnomADg,vcf,exact,0,AF,AF_NFE --check_existing --af --af_1kg --af_esp  --max_af --af_gnomad --fasta /media/cache/homo_sapiens/94_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz --variant_class --vcf_info_field ANN --plugin Gwava,region,/media/custom_data/gwava_scores/gwava_scores.bed.gz --plugin dbscSNV,/media/custom_data/dbscSNV/hg19/dbscSNV1.1_hg19.txt.gz --plugin ExACpLI,/media/custom_data/ExACpLI/ExACpLI_values.txt  --plugin GeneSplicer,/media/custom_data/genesplicer/bin/linux/genesplicer,/media/custom_data/genesplicer/human --plugin MaxEntScan,/media/custom_data/MaxEntScan --plugin SpliceRegion --plugin CADD,/media/custom_data/cadd/whole_genome_SNVs.tsv.gz,/media/custom_data/cadd/InDels.tsv.gz --plugin Condel,/media/custom_data/condel,b --plugin Carol --plugin dbNSFP,/media/custom_data/dbNSFP/hg19/dbNSFP_hg19.gz,GERP++_RS,phyloP20way_mammalian,SiPhy_29way_logOdds,phastCons20way_mammalian,MutationAssessor_score,SIFT_score,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,MutationTaster_score,MutationTaster_pred,MutationTaster_model,MutationTaster_AAE,FATHMM_score,FATHMM_pred,fathmm-MKL_coding_score,fathmm-MKL_coding_pred,PROVEAN_score,PROVEAN_pred,CADD_phred,DANN_score,Eigen-phred,Eigen-PC-phred,MetaSVM_score,MetaSVM_pred,MetaLR_score,MetaLR_pred,REVEL_score,GM12878_fitCons_score,H1-hESC_fitCons_score,HUVEC_fitCons_score_rankscore,LRT_score,LRT_pred,MutPred_score,VEST3_score,M-CAP_score --custom /media/custom_data/gerp/All_hg19_RS.bw,GERP,bigwig --custom /media/custom_data/phastcons100/hg19.100way.phastCons.bw,phastCons,bigwig --custom /media/custom_data/phyloP100/hg19.100way.phyloP100way.bw,phyloP,bigwig

#srun shifter --image=mcfonsecalab/variantutils:0.4 bcftools sort $output_vcf -Oz -o ${output_vcf/vcf.gz/sort.vcf.gz}
#zcat $output_vcf | grep "^#" > h && zcat $output_vcf | grep -v "^#" | sort -V -k1,1 -k2,2 | cat h - | shifter --image=ummidock/ubuntu_base:latest bgzip  > ${output_vcf/vcf/sort.vcf}
srun shifter --image=ummidock/ubuntu_base:latest tabix -p vcf $output_vcf
#srun shifter --image=ummidock/ubuntu_base:latest tabix -p vcf ${output_vcf/vcf.gz/sort.vcf.gz}

srun shifter --image=mcfonsecalab/variantutils:0.4 python ~/git_repos/bioinfo_utils/python-scripts/vcf-tools/prediction_tools/addReMM_scores.py $output_vcf ~/resources/predictionTools/reMM_genomiser/ReMM.v0.3.1.tsv.gz ${output_vcf/vcf.gz/ReMM.vcf} 

srun /home/pedro.barbosa/software/annovar/table_annovar.pl ${output_vcf/vcf.gz/ReMM.vcf} /home/pedro.barbosa/software/annovar/humandb/ -buildver hg19 -out ${output_vcf/vcf.gz/spidex} -remove -protocol spidex -operation f -nastring . -vcfinput

mv ${output_vcf/vcf.gz/spidex}*vcf ${output_vcf/.gz/}
srun shifter --image=ummidock/ubuntu_base:latest bgzip --force ${output_vcf/.gz/}
srun shifter --image=ummidock/ubuntu_base:latest tabix -p vcf ${output_vcf}
mv ${output_vcf}* $outdir
#cd ../ && rm -rf $SLURM_JOB_ID 

