#!/bin/bash
display_usage(){
    printf "
1st argument is the VCF file to annotate.
2nd argument is the output file. 
3rd argument is the output directory.
4th argument is optional. Refer to the tools to annotate, comma separated. Default: annotate all tools.
5th argument is optional. Create and run a slurm sbatch script on the fly. Default: true. Values: [true|false|-]. Set '-' to skip the argument.\n

Included tools:
gerp
phastcons
phyloP
fitcons
linsight
siphy
gwava
eigen
fathmmMKL
funseq
traP
spidex
traP
spliceai
"
}

if [[ -z $1 || -z $2 || -z $3 ]]; then
    printf "Error. Please set the required arguments for the script.\n"
    display_usage
    exit 1
fi
invcf=$(readlink -f $1)
outfile=$(readlink -f $2)
outdir=$(readlink -f $3)

annotate(){
    if [[ $1 == "spliceai" ]];then
#cat <<EOM >$PWD/custom.lua
#function round(num, numDecimalPlaces)
#  local mult = 10^(numDecimalPlaces or 0)
#  return math.floor(num * mult + 0.5) / mult
#end
#EOM
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/spliceai/whole_genome_filtered_spliceai_scores.vcf.gz"
fields=["SYMBOL","DS_AG","DS_AL","DS_DG","DS_DL","DP_AG","DP_AL","DP_DG","DP_DL"]
names=["gene_name","DS_AG","DS_AL","DS_DG","DS_DL","DP_AG","DP_AL","DP_DG","DP_DL"]
ops=["self","self","self","self","self","self","self","self","self"]

#[[postannotation]]
#name="SpliceAI"
#fields=["gene_name","DS_AG","DS_AL","DS_DG","DS_DL","DP_AG","DP_AL","DP_DG","DP_DL"]
#op="lua:join({gene_name, DS_AS, DS_AL, DS_DG, DS_DL, DP_AG, DP_AL, DP_DG, DP_DL}, '|')"
#type="String"
Description="Format: SpliceAIv1.2.1 variant annotation. These include delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL"

EOM
    elif [[ $1 == "gerp" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/gerp/gerp_hg19.bed.gz"
columns=[4]
names=["GERP"]
ops=["self"]

EOM
    
    elif [[ $1 == "phastcons" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]] 
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/phastcons100/phastCons_hg19.bed.gz"
columns=[4]
names=["phastCons"]
ops=["self"]

EOM

    elif [[ $1 == "phyloP" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/phyloP100/phyloP_hg19.bed.gz"
columns=[4]
names=["phyloP"]
ops=["self"]

EOM

   elif [[ $1 == "fitcons" ]]; then #nochr 
cat <<EOM >>$PWD/anno.conf 
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/wgsa_dbNSFP/resources/fitConsv1.01/fitcons_v1.01.bed.gz"
columns = [4]
names=["fitcons"]
ops=["self"]

EOM
    elif [[ $1 == "linsight" ]]; then #nochr
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/wgsa_dbNSFP/resources/LINSIGHT/LINSIGHT_hg19.bed.bgz"
columns=[4]
names=["linsight_g"]
ops=["self"]

EOM

     elif [[ $1 == "siphy" ]];then #nochr
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/wgsa_dbNSFP/resources/SiPhy/SiPhy_hg19.bed.gz"
columns=[5]
names=["SiPhy"]
ops=["self"]

EOM
    elif [[ $1 == "gwava" ]]; then #nochr
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/gwava_scores/gwava_scores.bed.gz"
columns=[5]
names=["Gwava"]
ops=["uniq"]

EOM
      elif [[ $1 == "ReMM" ]]; then #nochr #notworking
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/reMM_genomiser/ReMM.v0.3.1.txt.bgz"
columns=[3]
names=["ReMM"]
ops=["self"]

EOM
    elif [[ $1 == "eigen" ]];then #nochr
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/wgsa_dbNSFP/resources/Eigen/eigen_hg19.txt.gz"
columns=[7,9]
names=["Eigen","Eigen-PC"]
ops=["self","self"]

EOM
    elif [[ $1 == "DANN" ]];then #nochr noheader, lets see if there is problems. same problem as fathmmMKL
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/wgsa_dbNSFP/resources/DANN/DANN_hg19.txt.gz"
columns=[5]
names=["DANN"]
ops=["self"]

EOM
    elif [[ $1 == "fathmmMKL" ]];then #nochr
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/wgsa_dbNSFP/resources/fathmmMKL/fathmm-MKL_Current.txt.gz"
columns=[8]
names=["fathmmMKL"]
ops=["self"]

EOM

    elif [[ $1 == "funseq" ]]; then #nochr 
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/wgsa_dbNSFP/resources/funseq2/funseq216_hg19.bed.gz"
columns=[6]
names=["funseq2"]
ops=["self"]

EOM
    elif [[ $1 == "traP" ]]; then #nochr
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/traP/traP_v2.txt.gz"
columns=[6]
names=["traP"]
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
   fi
}

#ReMM
#spliceai
#DANN
#tools=(spliceai)
tools=(gerp phastcons phyloP fitcons linsight siphy gwava fathmmMKL eigen funseq traP spidex spliceai)
if [[ -f "$PWD/anno.conf" ]]; then
    rm "$PWD/anno.conf"
fi


if [[ $4 == "-" || -z "$4" ]]; then
    for elem in "${tools[@]}"
        do
        annotate $elem
        done
else
    IFS=','
    read -r -a array <<< "$4"
    for elem in "${array[@]}"
        do
            if [[ ! " ${tools[@]} " =~ " ${elem} " ]]; then
                printf "${elem} is not included in the list of available tools"
                display_usage
                exit 1
            else
                annotate $elem
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
    cmd="vcfanno -lua $lua_file $config $invcf | shifter --image=ummidock/ubuntu_base:latest bgzip  > $(basename $outfile)"
else
    cmd="vcfanno $config $invcf | shifter --image=ummidock/ubuntu_base:latest bgzip > $(basename $outfile)"
fi

   cat > $PWD/vcfAnno.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=vcfAnnot
#SBATCH --time=72:00:00
#SBATCH --mem=150G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
##SBATCH --workdir=/home/pedro.barbosa/scratch/vep
#SBATCH --output=%j_vcfanno.log
#SBATCH --image=mcfonsecalab/variantutils:0.5

cd /home/pedro.barbosa/scratch/vep
mkdir \${SLURM_JOB_ID}_vcfAnno && cd \${SLURM_JOB_ID}_vcfAnno
srun shifter $cmd
mv $(basename $outfile) $outdir
#cd ../ && rm -rf \$SLURM_JOB_ID
EOL

if [[ -z "$5" || $5 == "true" || $5 == "-" ]]; then 
    sbatch $PWD/vcfAnno.sbatch
fi


