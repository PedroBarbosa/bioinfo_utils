#!/bin/bash
display_usage(){
    printf "
1st argument is the VCF file to annotate.
2nd argument is the output file. 
3rd argument is the output directory.
4th argument is the genome version build of the input VCF. Default: hg19. Values:[hg19|hg38|-]
5th argument is optional. Refer to the tools to annotate, comma separated. Default: annotate all tools.\n

Included tools:
MVP (hg19)
primateAI (hg19,hg38)
cardioboost (hg19)
gerp (hg19, hg38)
phastcons (hg19, hg38)
phyloP (hg19, hg38)
fitcons (hg19)
linsight (hg19)
siphy (hg19)
cdts (hg19, hg38)
gwava (hg19)
eigen (hg19)
fathmmMKL (hg19)
remm (hg19)
dann (hg19)
funseq (hg19, hg38)
traP (hg19)
spidex (hg19)
scap (hg19)
spliceai (hg19, hg38)
gnomad_genomes (hg19, hg38)\n"
}

if [[ -z $1 || -z $2 || -z $3 ]]; then
    printf "Error. Please set the required arguments for the script.\n"
    display_usage
    exit 1
fi
invcf=$(readlink -f $1)
outfile=$(readlink -f $2)
outdir=$(readlink -f $3)
if [[ -z "$4" || "$4" == "-" || "$4" == "hg19" ]]; then
    genome_version="hg19"
elif [[ "$4" == "hg38" ]]; then
    genome_version="hg38"
else
    printf "Please set a valid value for the 4th argument.\n"
    display_usage
    exit 1
fi
  
annotate_hg19(){
    if [[ $1 == "spliceai" ]];then
#cat <<EOM >$PWD/custom.lua
#function round(num, numDecimalPlaces)
#  local mult = 10^(numDecimalPlaces or 0)
#  return math.floor(num * mult + 0.5) / mult
#end
#EOM
cat <<EOM >>$PWD/anno.conf
[[annotation]]
#file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/spliceai/whole_genome_filtered_spliceai_scores.vcf.gz"
file="mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/spliceai/v1.3/hg19/spliceai_scores.raw.snv.hg19.vcf.gz"
fields=["SpliceAI"]
names=["SpliceAI"]
ops=["self"]

#[[postannotation]]
#name="SpliceAI"
#fields=["gene_name","DS_AG","DS_AL","DS_DG","DS_DL","DP_AG","DP_AL","DP_DG","DP_DL"]
#op="lua:join({gene_name, DS_AS, DS_AL, DS_DG, DS_DL, DP_AG, DP_AL, DP_DG, DP_DL}, '|')"
#type="String"
Description="Format: SpliceAIv1.3 variant annotation. These include delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). Format: ALLELE|SYMBOL|STRAND|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL"

[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/spliceai/v1.3/hg19/spliceai_scores.raw.indel.hg19.vcf.gz"
fields=["SpliceAI"]
names=["SpliceAI_ind"]
ops=["self"]

EOM

    elif [[ $1 == "cdts" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/CDTS/hg19/CDTS_hg19_N11257all.bed.gz"
columns=[4,5]
names=["CDTS_score", "CDTS_percentile"]
ops=["self", "self"]

EOM

    elif [[ $1 == "gerp" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/gerp/hg19/gerp_hg19.bed.gz"
columns=[4]
names=["GERP"]
ops=["mean"]

EOM
    
    elif [[ $1 == "phastcons" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]] 
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/phastcons100/hg19/phastCons_hg19.bed.gz"
columns=[4]
names=["phastCons"]
ops=["mean"]

EOM

    elif [[ $1 == "phyloP" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/phyloP100/hg19/phyloP_hg19.bed.gz"
columns=[4]
names=["phyloP"]
ops=["mean"]

EOM

   elif [[ $1 == "fitcons" ]]; then #nochr 
cat <<EOM >>$PWD/anno.conf 
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/wgsa_dbNSFP/resources/fitConsv1.01/fitcons_v1.01.bed.gz"
columns = [4]
names=["fitcons"]
ops=["mean"]

EOM
    elif [[ $1 == "linsight" ]]; then #nochr
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/wgsa_dbNSFP/resources/LINSIGHT/LINSIGHT_hg19.bed.bgz"
columns=[4]
names=["linsight_g"]
ops=["mean"]

EOM

     elif [[ $1 == "siphy" ]];then #nochr
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/wgsa_dbNSFP/resources/SiPhy/SiPhy_hg19.bed.gz"
columns=[5]
names=["SiPhy"]
ops=["mean"]

EOM
    elif [[ $1 == "gwava" ]]; then #nochr
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/gwava_scores/gwava_scores.bed.gz"
columns=[5,5]
names=["GWAVA","GWAVA_mean"]
ops=["uniq","mean"]

EOM
      elif [[ $1 == "remm" ]]; then #nochr #notworking
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/reMM_genomiser/ReMM.v0.3.1.txt.bgz"
columns=[3]
names=["ReMM"]
ops=["mean"]

EOM
    elif [[ $1 == "eigen" ]];then #nochr
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/wgsa_dbNSFP/resources/Eigen/eigen_hg19.txt.gz"
columns=[7,9]
names=["Eigen","Eigen-PC"]
ops=["self","self"]

EOM
    elif [[ $1 == "dann" ]];then #nochr noheader, lets see if there is problems. same problem as fathmmMKL
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/wgsa_dbNSFP/resources/DANN/DANN_hg19.txt.gz"
columns=[5]
names=["DANN"]
ops=["self"]

EOM

    elif [[ $1 == "fathmmMKL" ]]; then #nochr 
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/wgsa_dbNSFP/resources/fathmmMKL/hg19/fathmm-MKL_Current.txt.gz"
columns=[6]
names=["fathmmMKL"]
ops=["self"]

EOM

    elif [[ $1 == "funseq" ]]; then #nochr 
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/wgsa_dbNSFP/resources/funseq2/hg19/funseq216_hg19.bed.gz"
columns=[6,6]
names=["funseq2","funseq2_mean"]
ops=["self","mean"]

EOM

    elif [[ $1 == "traP" ]]; then #nochr
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/traP/v3/TraP_v3.txt.gz"
columns=[6]
names=["TraP"]
ops=["self"]

EOM
    elif [[ $1 == "spidex" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/spidex/spidex_hg19.txt.gz"
columns=[6,7]
names=["dpsi_max_tissue","dpsi_zscore"]
ops=["self","self"]

EOM

    elif [[ $1 == "scap" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/S-CAP/scap_COMBINED_v1.0.vcf.gz"
fields=["SCAP","SCAP"]
names=["SCAP","SCAP"]
ops=["self","self"]

EOM

   elif [[ $1 == "MVP" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/MVP/hg19/MVP_hg19.tsv.gz"
columns=[9]
names=["MVP"]
ops=["self"]

EOM

    elif [[ $1 == "primateAI" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/primateAI/hg19/primateAI_scores_v0.2.tsv.gz"
columns=[11]
names=["PrimateAI"]
ops=["self"]

EOM

   elif [[ $1 == "cardioboost" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/CardioBoost/hg19/CardioBoost_hg19.tsv.gz"
columns=[8]
names=["CardioBoost"]
ops=["self"]

EOM

   elif [[ $1 == "gnomad_genomes" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/gnomAD/hg19/gnomAD_v2.1_justImportantFields.vcf.gz"
fields=["AF","AF_nfe","AF_nfe_nwe"]
names=["gnomADg_AF", "gnomADg_AF_nfe", "gnomADg_AF_nwe"]
ops=["self", "self", "self"]

EOM
   fi
}


annotate_hg38(){
    if [[ $1 == "spliceai" ]];then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/spliceai/v1.3/hg38/spliceai_scores.raw.snv.hg38.vcf.gz"
fields=["SpliceAI"]
names=["SpliceAI"]
ops=["self"]

[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/spliceai/v1.3/hg38/spliceai_scores.raw.indel.hg38.vcf.gz"
fields=["SpliceAI"]
names=["SpliceAI_ind"]
ops=["self"]

EOM

    elif [[ $1 == "cdts" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/CDTS/hg38/CDTS_hg38_N11257all.bed.gz"
columns=[4,5]
names=["CDTS_score", "CDTS_percentile"]
ops=["self", "self"]

EOM


    elif [[ $1 == "gerp" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/gerp/hg38/gerp_hg38.bed.gz"
columns=[4]
names=["GERP"]
ops=["mean"]

EOM

    elif [[ $1 == "phastcons" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]] 
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/phastcons100/hg38/phastCons_hg38.bed.gz"
columns=[4]
names=["phastCons"]
ops=["mean"]

EOM

    elif [[ $1 == "phyloP" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/phyloP100/hg38/phyloP_hg38.bed.gz"
columns=[4]
names=["phyloP"]
ops=["mean"]

EOM

    elif [[ $1 == "primateAI" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/primateAI/hg38/primateAI_scores_v0.2.tsv.gz"
columns=[11]
names=["PrimateAI"]
ops=["self"]

EOM

    elif [[ $1 == "gnomad_genomes" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
#with chrom
#file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/gnomAD/hg38/gnomAD_genomes_v3.0_with_chr.vcf.bgz"
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/gnomAD/hg38/gnomAD_genomes_v3.0_no_chr.vcf.bgz"
fields=["AF","AF_nfe"]
names=["gnomADg_AF", "gnomADg_AF_nfe"]
ops=["self", "self"]

EOM
   fi
}

tools=(gerp phastcons phyloP fitcons cdts linsight siphy gwava fathmmMKL eigen remm dann funseq traP spidex spliceai scap MVP cardioboost primateAI gnomad_genomes)
do_spliceAI="false"

if [[ -f "$PWD/anno.conf" ]]; then
    rm "$PWD/anno.conf"
fi


if [[ $5 == "-" || -z "$5" ]]; then
    do_spliceAI="true"
    for elem in "${tools[@]}"
        do
        if [[ $genome_version == "hg19" ]]; then
            annotate_hg19 $elem
        elif [[ $genome_version == "hg38" ]]; then
            annotate_hg38 $elem
        fi
        done
else
    IFS=','
    read -r -a array <<< "$5"
    for elem in "${array[@]}"
        do
            if [[ ! " ${tools[@]} " =~ " ${elem} " ]]; then
                printf "${elem} is not included in the list of available tools"
                display_usage
                exit 1
            else
                if [[ "${elem}" == "spliceai" ]]; then
                    do_spliceAI="true"
                fi

                if [[ $genome_version == "hg19" ]]; then
                    annotate_hg19 $elem
                elif [[ $genome_version == "hg38" ]]; then
                    annotate_hg38 $elem
                fi
                
	    fi
        done
fi

config=$(readlink -f $PWD/anno.conf)
if [ -f $PWD/custom.lua ]; then
    custom_lua=True
    lua_file="$PWD/custom.lua"
else
    custom_lua=False
fi

if [ $custom_lua == "True" ]; then
    cmd="vcfanno -lua $lua_file $config $invcf | shifter --image=ummidock/ubuntu_base:latest bgzip > ${outfile}_tmp"
else
    cmd="vcfanno $config $invcf | shifter --image=ummidock/ubuntu_base:latest bgzip > ${outfile}_tmp"
fi

cat > $PWD/vcfAnno.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=vcfAnnot
#SBATCH --time=72:00:00
#SBATCH --mem=120G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --output=%j_vcfanno.log
#SBATCH --image=mcfonsecalab/variantutils:latest 

srun shifter $cmd

if [[ $do_spliceAI == "true" ]]; then
    srun shifter --image=ummidock/ubuntu_base:latest tabix --force -p vcf ${outfile}_tmp
    srun zcat ${outfile}_tmp | shifter python /home/pedro.barbosa/git_repos/bioinfo_utils/python-scripts/vcf-tools/prediction_tools/split_SpliceAI_field.py - | shifter --image=ummidock/ubuntu_base:latest bgzip > $outfile
    srun rm ${outfile}_tmp*

else
    mv ${outfile}_tmp ${outfile}
fi

srun shifter --image=ummidock/ubuntu_base:latest tabix --force -p vcf $outfile
EOL

sbatch $PWD/vcfAnno.sbatch
