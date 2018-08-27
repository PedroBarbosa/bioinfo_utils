#!/bin/bash
#SBATCH --job-name=rs2VCF
#SBATCH --time=72:00:00
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --workdir=/home/pedro.barbosa/scratch/vep
#SBATCH --output=/home/pedro.barbosa/scratch/vep/%j_vep_cache.log

mkdir $SLURM_JOB_ID && cd $SLURM_JOB_ID
rsIDs=$(readlink -f $1) 
output_vcf=$(readlink -f $2)
outdir=$(readlink -f $3)

srun shifter -V=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.barbosa/resources/ensembl_vep:/media --image=ensemblorg/ensembl-vep:latest vep -i $rsIDs -o $output_vcf --format id  --hgvs --hgvsg --pick --compress_output bgzip  --cache --dir /media/cache --vcf --numbers --regulatory --force_overwrite --port 3337 -a GRCh37 --variant_class --vcf_info_field ANN --plugin Gwava,region,/media/custom_data/gwava_scores/gwava_scores.bed.gz --plugin dbscSNV,/home/pedro.barbosa/resources/ensembl_vep/custom_data/dbscSNV/hg19/dbscSNV1.1_hg19.txt.gz --plugin ExACpLI,/media/custom_data/ExACpLI/ExACpLI_values.txt --plugin GeneSplicer,/media/custom_data/genesplicer/bin/linux/genesplicer,/media/custom_data/genesplicer/human --plugin MaxEntScan,/media/custom_data/MaxEntScan --plugin SpliceRegion --plugin Condel,/home/pedro.barbosa/resources/ensembl_vep/custom_data/condel,b --plugin Carol --plugin CADD,/media/custom_data/cadd/whole_genome_SNVs.tsv.gz,/media/custom_data/cadd/InDels.tsv.gz --plugin dbNSFP,/media/custom_data/dbNSFP/hg19/dbNSFP_hg19.gz,GERP++_RS,phyloP20way_mammalian,SiPhy_29way_logOdds,phastCons20way_mammalian,MutationAssessor_score,SIFT_score,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,MutationTaster_score,MutationTaster_pred,MutationTaster_model,MutationTaster_AAE,FATHMM_score,FATHMM_pred,fathmm-MKL_coding_score,fathmm-MKL_coding_pred,PROVEAN_score,PROVEAN_pred,CADD_phred,DANN_score,Eigen-phred,Eigen-PC-phred,MetaSVM_score,MetaSVM_pred,MetaLR_score,MetaLR_pred,REVEL_score,GM12878_fitCons_score,H1-hESC_fitCons_score,HUVEC_fitCons_score_rankscore,LRT_score,LRT_pred,MutPred_score,VEST3_score,M-CAP_score

zcat $output_vcf | grep "^#" > h && zcat $output_vcf | grep -v "^#" | sort -V -k1,1 -k2,2 | cat h - | shifter --image=ummidock/ubuntu_base:latest bgzip  > ${output_vcf/vcf/sort.vcf}
shifter --image=ummidock/ubuntu_base:latest tabix -p vcf ${output_vcf/vcf/sort.vcf}

srun shifter --image=mcfonsecalab/variantutils:0.4 python ~/git_repos/bioinfo_utils/python-scripts/vcf-tools/prediction_tools/addReMM_scores.py ${output_vcf/vcf/sort.vcf} ~/resources/predictionTools/reMM_genomiser/ReMM.v0.3.1.tsv.gz ${output_vcf/vcf.gz/ReMM.vcf} 

srun /home/pedro.barbosa/software/annovar/table_annovar.pl ${output_vcf/vcf.gz/ReMM.vcf} /home/pedro.barbosa/software/annovar/humandb/ -buildver hg19 -out ${output_vcf/vcf.gz/spidex} -remove -protocol spidex -operation f -nastring . -vcfinput

mv * ../${SLURM_JOB_ID}* $outdir
