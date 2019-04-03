#!/bin/bash
set -e
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
           11th argument is optional. Refers whether we should parallelize VEP run per chromosome. Values:[true|false|-]. Set '-' to skip the argument.
           12th argument is optional. Refers to any additional argument that will compose the VEP command.
\n" 
           
}

if [[ -z "$1" || -z "$2" || -z "$3" || -z "$4"  || -z "$5" ]] ; then
    printf "ERROR.Please set the required arguments.\n"
    display_usage
    exit 1
fi

BCFTOOLS="shifter --image=mcfonsecalab/variantutils:0.5 bcftools"
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
if [[ ! -d "$OUT_DIR" ]]; then
    mkdir $OUT_DIR
fi
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
if [ -z "$7" -o "$7" == "-" -o  "$7" == "true" ] && [ "$genome_version" != "mm10" ]; then
    #MISSING CUSTOM SCORES FOR HG38"
    if [[ "$genome_version" == "hg38" ]]; then
        printf "INFO. Custom plugins can't be used as prediction tools are not ready (yet) for the latest version. Skipping this flag.\n"
    else
        runTools="true"
        BASE_CMD="$BASE_CMD --plugin dbscSNV,/media/custom_data/dbscSNV/hg19/dbscSNV1.1_hg19.txt.gz \
--plugin ExACpLI,/media/custom_data/ExACpLI/ExACpLI_values.txt \
--plugin GeneSplicer,/media/custom_data/genesplicer/bin/linux/genesplicer,/media/custom_data/genesplicer/human \
--plugin MaxEntScan,/media/custom_data/MaxEntScan \
--plugin SpliceRegion \
--plugin Condel,/media/custom_data/condel,b \
--plugin Carol \
--plugin CADD,/media/custom_data/cadd/hg19/v1.4/whole_genome_SNVs.tsv.gz,/media/custom_data/cadd/hg19/v1.4/InDels.tsv.gz \
--plugin dbNSFP,'consequence=ALL',/media/custom_data/dbNSFP/hg19/dbNSFP4.0b1a_hg19.txt.gz,GERP++_RS,phyloP30way_mammalian,29way_logOdds,phastCons30way_mammalian,MutationAssessor_score,SIFT4G_score,SIFT_score,Polyphen2_HDIV_score,Polyphen2_HVAR_score,MutationTaster_score,FATHMM_score,fathmm-MKL_coding_score,PROVEAN_score,CADD_phred,DANN_score,Eigen-pred_coding,Eigen-PC-phred_coding,MetaSVM_score,MetaSVM_pred,MetaLR_score,MetaLR_pred,REVEL_score,integrated_fitCons_score,LRT_score,LRT_pred,MutPred_score,VEST4_score,M-CAP_score,LINSIGHT \
--custom /media/custom_data/gerp/All_hg19_RS.bw,GERP_vep,bigwig,exact,0 \
--custom /media/custom_data/gerp/All_hg19_RS.bw,GERP_vep_all,bigwig,overlap,0 \
--custom /media/custom_data/phyloP100/hg19.100way.phyloP100way.bw,phyloP_vep,bigwig,exact,0 \
--custom /media/custom_data/phastcons100/hg19.100way.phastCons.bw,phastCons_vep,bigwig,exact,0 \
--custom /media/custom_data/gnomAD/hg19/gnomAD_v2.1_justImportantFields.vcf.gz,gnomADg,vcf,exact,0,AF,AF_nfe,AF_nfe_nwe" # \
#--plugin Gwava,region,/media/custom_data/gwava_scores/gwava_scores.bed.gz \
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


#### ADITIONAL ARGS ####
if [[ ! -z "${12}" ]]; then
    printf "INFO. ${12} string will be passed to the main command.\n"
    BASE_CMD="$BASE_CMD ${12}"
fi

BASE_CMD="$BASE_CMD -i $IN_VCF"
PARALLEL=false
FINAL_OUT=$OUT
##### OUTPUT #####
if [[ "${11}" != "true" ]]; then
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

    cat > $WORKDIR/runVEP.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=runVEP
#SBATCH --time=72:00:00
#SBATCH --mem=150G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=%j_runvep.log

cd $WORKDIR && mkdir \${SLURM_JOB_ID}_noParallel && cd \${SLURM_JOB_ID}_noParallel
$BASE_CMD
srun shifter --image=ummidock/ubuntu_base:latest tabix -p vcf ${OUT}.vcf.bgz
mv ${OUT}.vcf.bgz* $OUT_DIR
#cd ../ && rm -rf \$SLURM_JOB_ID*
EOL

#### PARALLEL RUN ####
elif [[ "${11}" == "true" ]]; then
    printf "Extracting chromosome names to parallelize the run..\n"
    if [[ ${IN_VCF: -4} == ".vcf" || ${IN_VCF: -7} == "vcf.bgz" || ${IN_VCF: -6} == "vcf.gz" ]]; then

        $BCFTOOLS view -H $IN_VCF | awk '{print $1}' | sort -V | uniq > $WORKDIR/listChroms.txt 
        mapfile -t chroms < <( $BCFTOOLS view -H $IN_VCF | awk '{print $1}' | sort -V | uniq )    
    else
        awk '{print $1}' $IN_VCF | sort -V | uniq > $WORKDIR/listChroms.txt
        mapfile -t chroms < <( awk '{print $1}' $IN_VCF | sort -V | uniq )
    fi
    number_chroms="$((${#chroms[@]}-1))"
    if [[ "$number_chroms" -lt 15 ]]; then
        n_arrays=$number_chroms
    else
        n_arrays=15
    fi 
    PARALLEL=true
    OUT="\${chrom_list[\$SLURM_ARRAY_TASK_ID]}"

    cat << EOF > $WORKDIR/runVEP.sbatch
#!/bin/bash
#SBATCH --job-name=runVEP
#SBATCH --time=72:00:00
#SBATCH --mem=30G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=%j_runvep.log
#SBATCH --array=0-$number_chroms%$n_arrays

cd $WORKDIR 
readarray -t chrom_list < listChroms.txt
mkdir \${SLURM_JOB_ID}_withParallel && cd \${SLURM_JOB_ID}_withParallel
$BASE_CMD --chr \${chrom_list[\$SLURM_ARRAY_TASK_ID]} --output_file \${chrom_list[\$SLURM_ARRAY_TASK_ID]}.vcf.bgz --vcf --vcf_info_field ANN --compress_output bgzip
srun shifter --image=ummidock/ubuntu_base:latest tabix -p vcf \${chrom_list[\$SLURM_ARRAY_TASK_ID]}.vcf.bgz
mv \${chrom_list[\$SLURM_ARRAY_TASK_ID]}.vcf.bgz* ../
#cd ../ && rm -rf \$SLURM_JOB_ID*
EOF

elif [[ "${11}" != "false" && "$11" != "-" ]]; then
    printf "Please set a valid value for the 11th argument.\n"
    exit 1
fi

printf "\n$BASE_CMD\n"
job_submission_vep=$(sbatch $WORKDIR/runVEP.sbatch)
job_id_vep=${job_submission_vep##* }
if [[ $PARALLEL == "true" ]]; then
    echo -e "#!/bin/bash
#SBATCH --job-name=mergeVCFsv
#SBATCH --time=72:00:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=mcfonsecalab/variantutils:0.5
#SBATCH --workdir=$WORKDIR
#SBATCH --output=%j_concatvep.log\n" > $WORKDIR/concatVCFs_vep.sbatch
    echo -e "srun shifter bcftools concat -a -Oz *.vcf.bgz | shifter bcftools sort -Oz -o ${FINAL_OUT}.vcf.bgz\n" >> $WORKDIR/concatVCFs_vep.sbatch
    echo -e "srun shifter --image=ummidock/ubuntu_base:latest tabix -p vcf ${FINAL_OUT}.vcf.bgz\n" >> $WORKDIR/concatVCFs_vep.sbatch
    echo -e "mv ${FINAL_OUT}* $OUT_DIR" >> $WORKDIR/concatVCFs_vep.sbatch
    job_submission_concat=$(sbatch --depend=afterok:$job_id_vep $WORKDIR/concatVCFs_vep.sbatch)
    job_id_concat=${job_submission_concat##* }
fi

if [[ "$runTools" == "true" ]]; then
    cat << EOF >> $WORKDIR/runAdditionalTools.sbatch 
#!/bin/bash
#SBATCH --job-name=runAdditionalTools
#SBATCH --time=72:00:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=mcfonsecalab/variantutils:0.5
#SBATCH --output=%j_runAdditionalTools.log

cd $WORKDIR && mkdir \${SLURM_JOB_ID}_additionalTools && cd \${SLURM_JOB_ID}_additionalTools
if [[ -f "$OUT_DIR/${FINAL_OUT}.vcf.bgz" ]]; then
    if [[ $genome_version == "hg19" ]]; then
        printf "Adding ReMM and DANN scores..\n"
        srun shifter --image=mcfonsecalab/variantutils:0.5 python ~/git_repos/bioinfo_utils/python-scripts/vcf-tools/prediction_tools/add_prediction_scores.py $OUT_DIR/${FINAL_OUT}.vcf.bgz ${FINAL_OUT}_ReMM_DANN.vcf
        #/home/pedro.barbosa/software/miniconda3/bin/spliceai
        #srun cat ${FINAL_OUT}_ReMM_DANN.vcf | shifter --image=mcfonsecalab/variantutils:0.5 spliceai -R genome.fa -A grch37 | shifter --image=ummidock/ubuntu_base:latest bgzip > ${FINAL_OUT}_ReMM_DANN_spliceAI.vcf.gz
        srun shifter --image=ummidock/ubuntu_base:latest bgzip --force ${FINAL_OUT}_ReMM_DANN.vcf

        printf "Running vcfanno to add remaining scores..\n"
        run_vcfanno=\$(/home/pedro.barbosa/git_repos/bioinfo_utils/sbatch-scripts/variant_annotation/prediction_tools/vcfAnno.sh \$PWD/${FINAL_OUT}_ReMM_DANN.vcf.gz ${FINAL_OUT}_final.vcf.bgz $OUT_DIR - true)
        echo \$run_vcfanno
        \$run_vcfanno
        mv *sbatch ../*sbatch $OUT_DIR
# && cd ../ # && rm -rf \${SLURM_JOB_ID}_additionalTools

    else
        printf "Additional tools can't be run (yet) on the latest genome build. Skipping this step.\n"
    fi
else
    printf "$OUT_DIR/${FINAL_OUT}.vcf.bgz does not exist. Additional scores won't be added.\n"
fi
EOF

    if [[ "$PARALLEL" == "true" ]]; then
        sbatch --depend=afterok:$job_id_concat $WORKDIR/runAdditionalTools.sbatch
    else
        sbatch --depend=afterok:$job_id_vep $WORKDIR/runAdditionalTools.sbatch
    fi
 
fi
