#!/bin/bash
display_usage(){
   printf "
           1st argument is the VCF/rsIDs/hgvs/vep_default file to annotate.
	   2nd argument is the basename for the output file (vcf.gz extension will be automatically added).
           3rd argument is the output directory.
           4th argument is the genome version to use. Values:[hg19,hg38,mm10]. 
           5th argument is the cache version to use. Values:[95,94].
           6th argument is optional. Refers to the set of annotations to use. Default:ensembl. Values:[ensembl|refseq|merged|-]. Set '-' to skip the argument.
           7th argument is optional. Refers whether custom annotations (e.g. gnomAD genomes frequencies, dbNSFP, conservation scores) should be added. Default:true. Values:[true|false|-]. Set '-' to skip the argument.
           8th argument is optional. Refers whether allele frequencies should be added. Only work for human caches. Default:true. Values:[true|false|-]. Set '-' to skip the argument.
           9th argument is optional. Refers whether variant consequences should be picked. Argument should be set within quotes, so all the pick option are automatically passed to the VEP main command. Default:false. Values: [false|-|'command']. Set '-' to skip the argument. 'command' value may include any of the following flags:'--pick --pick_allele --per_gene --pick_order tsl,appris,rank'
           10th argument is optional. Refers to the output format. Default:vcf. Values:[vcf|tab|json|-].
           11th argument is optional. Refers to any additional argument that will compose the VEP command.
\n" 
           
}

if [[ -z "$1" || -z "$2" || -z "$3" || -z "$4"  || -z "$5" ]] ; then
    printf "ERROR.Please set the required arguments.\n"
    display_usage
    exit 1
fi

WORKDIR="/home/pedro.barbosa/scratch/vep"
if [[ ! -f $(readlink -f "$1") ]]; then
    printf "ERROR. Please set a valid file in the 1st argument.\n"
    display_usage
    exit 1
else
    IN_VCF=$(readlink -f "$1")
fi
OUT=$2
OUT_DIR=$(readlink -f "$3")
BASE_CMD="srun shifter -V=/mnt/nfs/lobo/IMM-NFS/ensembl_vep:/media --image=ensemblorg/ensembl-vep:latest vep \
--force_overwrite --stats_text --hgvs --hgvsg --pubmed --check_existing \
--cache --dir /media/cache --offline --sift b --polyphen b --numbers --regulatory --variant_class"

genomes=(hg19 hg38 mm10)
cache=(95 94)
annotations=(ensembl refseq merged)
genome_version="$4"
cache_version="$5"
if [[ -z "$6" || "$6" == "-" ]];then
    annotation="ensembl"
else
    annotation="$6"
fi

if [[ ! " ${annotations[@]} " =~ " $annotation " ]]; then
    printf "ERROR. $annotation is not included in the list of available annotations.\n"
    display_usage
    exit 1
elif [[ ! " ${genomes[@]} " =~ " $genome_version " ]]; then
    printf "ERROR. $genome_version is not included in the list of available genomes.\n"
    display_usage
    exit 1
elif [[ ! " ${cache[@]} " =~ " $cache_version " ]]; then
    printf "ERROR. $cache_version is not included in the list of available caches.\n"
    display_usage
    exit 1
fi

if [[ "$genome_version" == "mm10" ]]; then
    printf "INFO. When mus musculus cache is used, it does not matter the cache version and the set of annotations set. VEP defaults to the 95 cache, ensembl annotations.\n"
    ASSEMBLY="GRCm38"
    FASTA="/media/cache/mus_musculus/95_GRCm38/Mus_musculus.GRCm38.dna.toplevel.fa.gz"

elif [[ "$annotation" == "ensembl" ]]; then
    annot_dir="/media/cache/homo_sapiens"
elif [[ "$annotation" == "refseq" ]]; then
    BASE_CMD="$BASE_CMD --refseq"
    annot_dir="/media/cache/homo_sapiens_refseq"
else
    BASE_CMD="$BASE_CMD --merged"
    annot_dir="/media/cache/home_sapiens_merged"
fi

if [[ "$genome_version" == "hg19" ]]; then
    ASSEMBLY="GRCh37"
    if [[ "$cache_version" == "94" ]]; then
        FASTA="$annot_dir/94_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"
    elif [[ "$cache_version" == "95" ]]; then
        FASTA="$annot_dir/95_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"
    fi

elif [[ "$genome_version" == "hg38" ]];then
    ASSEMBLY="GRCh38"
    if [[ "$cache_version" == "94" ]]; then
        FASTA="$annot_dir/94_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
    elif [[ "$cache_version" == "95" ]]; then
        FASTA="$annot_dir/95_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
    fi
fi
BASE_CMD="$BASE_CMD --cache_version $cache_version -a $ASSEMBLY --fasta $FASTA"

###### PLUGINS #######
# [  "$A" -eq "0" -o "$B" -ne "0"  ] && [ "$C" -eq "0" ];
if [ -z "$7" -o "$7" == "-" -o  "$7" == "true" ] && [ "$genome_version" != "mm10" ]; then
    #MISSING CUSTOM SCORES FOR HG38"
    if [[ "$genome_version" == "hg38" ]]; then
        printf "INFO. Custom plugins can't be used as prediction tools are not ready (yet) for the latest version. Skipping this flag.\n"
    else
        runTools="true"
        BASE_CMD="$BASE_CMD --plugin Gwava,region,/media/custom_data/gwava_scores/gwava_scores.bed.gz \
--plugin dbscSNV,/media/custom_data/dbscSNV/hg19/dbscSNV1.1_hg19.txt.gz \
--plugin ExACpLI,/media/custom_data/ExACpLI/ExACpLI_values.txt \
--plugin GeneSplicer,/media/custom_data/genesplicer/bin/linux/genesplicer,/media/custom_data/genesplicer/human \
--plugin MaxEntScan,/media/custom_data/MaxEntScan \
--plugin SpliceRegion \
--plugin Condel,/media/custom_data/condel,b \
--plugin Carol \
--plugin CADD,/media/custom_data/cadd/hg19/v1.3/whole_genome_SNVs.tsv.gz,/media/custom_data/cadd/hg19/v1.3/InDels.tsv.gz \
--plugin dbNSFP,'consequence=ALL',/media/custom_data/dbNSFP/hg19/dbNSFP4.0b1a_hg19.txt.gz,GERP++_RS,phyloP30way_mammalian,29way_logOdds,phastCons30way_mammalian,MutationAssessor_score,SIFT4G_score,SIFT_score,Polyphen2_HDIV_score,Polyphen2_HVAR_score,MutationTaster_score,FATHMM_score,fathmm-MKL_coding_score,PROVEAN_score,CADD_phred,DANN_score,Eigen-pred_coding,Eigen-PC-phred_coding,MetaSVM_score,MetaSVM_pred,MetaLR_score,MetaLR_pred,REVEL_score,integrated_fitCons_score,LRT_score,LRT_pred,MutPred_score,VEST4_score,M-CAP_score,LINSIGHT \
--custom /media/custom_data/gnomeAD/hg19/gnomAD_v2.1_justImportantFields_chrRenamed.vcf.gz,gnomADg,vcf,exact,0,AF,AF_nfe,AF_nfe_nwe \
--custom /media/custom_data/gerp/All_hg19_RS.bw,GERP,bigwig \
--custom /media/custom_data/phastcons100/hg19.100way.phastCons.bw,phastCons,bigwig \
--custom /media/custom_data/phyloP100/hg19.100way.phyloP100way.bw,phyloP,bigwig \
--custom /media/custom_data/wgsa_dbNSFP/resources/LINSIGHT/LINSIGHT_hg19.bw,LINSIGHTg,bigwig"
    fi
elif [ "$genome_version" == "mm10" ]; then
    printf "INFO. Plugins can't be set when running for mus musculus cache. Skipping this flag whatever the value was.\n"

elif [ "$7" != "false" ]; then
    printf "ERROR. Please set a valid value for the 7th argument.\n"
    display_usage
    exit 1
else
    printf "INFO. Plugins were not set.\n"
fi


##### ALLELE FREQUENCIES #####
if [ -z "$8" -o "$8" == "-" -o "$8" == "true" ] && [  "$genome_version" != "mm10" ]; then
    BASE_CMD="$BASE_CMD --af --af_1kg --af_esp  --max_af --af_gnomad"
elif [ "$genome_version" == "mm10" ]; then
    printf "INFO. Frequencies can't be set when running with mus musculus cache. Skipping this flag.\n"
elif [ "$8" != "false" ]; then
    printf "ERROR. Please set a valid value for the 8th argument.\n"
    display_usage
    exit 1
else
    printf "INFO. Frequencies were not set.\n"
fi

##### PICK OPTIONS #####
if [[ ! -z "$9" && "$9" != "-" ]]; then
    printf "INFO. Pick options set.\n"
    BASE_CMD="$BASE_CMD $9"
fi

##### OUTPUT #####
if [[ -z "${10}" || ${10} == "vcf" || ${10} == "-" ]]; then
    BASE_CMD="$BASE_CMD --vcf --vcf_info_field ANN --compress_output bgzip --output_file ${OUT}.vcf.bgz"
elif [[ ${10} == "tab" ]]; then
    BASE_CMD="$BASE_CMD --tab --compress_output gzip --output_file ${OUT}.tsv.gz"
elif [[ ${10} == "json" ]]; then
    BASE_CMD="$BASE_CMD --json --output_file ${OUT}.json"
else
    printf "ERROR. Please set a valid output format in the 10th argument.\n"
    display_usage
    exit 1
fi

#### ADITIONAL ARGS ####
if [[ ! -z "${11}" ]]; then
    printf "INFO. ${11} string will be passed to the main command.\n"
    BASE_CMD="$BASE_CMD ${11}"
fi

BASE_CMD="$BASE_CMD -i $IN_VCF"
printf "\n$BASE_CMD\n"
cat > $WORKDIR/runVEP.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=runVEP
#SBATCH --time=72:00:00
#SBATCH --mem=150G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=%j_vep.log

cd $WORKDIR && mkdir \$SLURM_JOB_ID && cd \$SLURM_JOB_ID

$BASE_CMD
#srun shifter --image=mcfonsecalab/variantutils:0.4 bcftools sort $output_vcf -Oz -o ${output_vcf/vcf.gz/sort.vcf.gz}
#zcat $output_vcf | grep "^#" > h && zcat $output_vcf | grep -v "^#" | sort -V -k1,1 -k2,2 | cat h - | shifter --image=ummidock/ubuntu_base:latest bgzip  > ${output_vcf/vcf/sort.vcf}
#srun shifter --image=ummidock/ubuntu_base:latest tabix -p vcf $output_vcf
#srun shifter --image=ummidock/ubuntu_base:latest tabix -p vcf ${output_vcf/vcf.gz/sort.vcf.gz}
if [[ -f "${OUT}.vcf.bgz" && "$runTools" == "true" ]]; then
    vep_out="${OUT}.vcf.bgz"
    srun shifter --image=ummidock/ubuntu_base:latest tabix -p vcf \$vep_out 
    srun shifter --image=mcfonsecalab/variantutils:0.4 python ~/git_repos/bioinfo_utils/python-scripts/vcf-tools/prediction_tools/addReMM_scores.py ${OUT}.vcf.bgz /mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/reMM_genomiser/ReMM.v0.3.1.txt.gz \${vep_out/vcf.bgz/ReMM.vcf} 
    srun /home/pedro.barbosa/software/annovar/table_annovar.pl \${vep_out/vcf.bgz/ReMM.vcf} /home/pedro.barbosa/software/annovar/humandb/ -buildver hg19 -out \${vep_out/vcf.bgz/spidex} -remove -protocol spidex -operation f -nastring . -vcfinput
    mv \${vep_out/vcf.bgz/spidex}*vcf \${vep_out/.bgz/}
    srun shifter --image=mcfonsecalab/variantutils:latest spliceai -I \${vep_out/.bgz/} -R /mnt/nfs/lobo/IMM-NFS/genomes/hg19/Sequence/WholeGenomeFasta/genome.fa -A grch37 | \
shifter --image=ummidock/ubuntu_base:latest bgzip --force > \${vep_out/.bgz/.gz}
    srun shifter --image=ummidock/ubuntu_base:latest tabix -p vcf \${vep_out/.bgz/.gz}
    mv \${vep_out/.bgz/.gz} $OUT_DIR
else
    mv ${OUT}* $OUT_DIR
fi
cd ../ && rm -rf \$SLURM_JOB_ID 
EOL
sbatch $WORKDIR/runVEP.sbatch
