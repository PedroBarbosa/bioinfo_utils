#!/bin/bash
display_usage(){
    printf "
1st argument is the VCF file to annotate.
2nd argument is the output file. 
3rd argument is the output directory.
4th argument is optional. Refer to the tools to annotate, comma separated. Default: annotate all tools.\n
Included tools:
fitcons
linsight
ReMM
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
cat <<EOM >$PWD/custom.lua
function round(num, numDecimalPlaces)
  local mult = 10^(numDecimalPlaces or 0)
  return math.floor(num * mult + 0.5) / mult
end
EOM
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/home/pedro.barbosa/spliceai/vcfanno/whole_genome_filtered_spliceai_scores_withChr.vcf.gz"
fields=["SYMBOL","DS_AG","DS_AL","DS_DG","DS_DL","DP_AG","DP_AL","DP_DG","DP_DL"]
names=["gene_name","DS_AG","DS_AL","DS_DG","DS_DL","DP_AG","DP_AL","DP_DG","DP_DL"]
ops=["self","self","self","self","self","self","self","self","self"]

[[postannotation]]
name="SpliceAI"
fields=["gene_name","DS_AG","DS_AL","DS_DG","DS_DL","DP_AG","DP_AL","DP_DG","DP_DL"]
op="lua:join({alt[1],gene_name, DS_AS, DS_AL, DS_DG, DS_DL, DP_AG, DP_AL, DP_DG, DP_DL}, '|')"
type="String"
Description="Format: SpliceAIv1.2.1 variant annotation. These include delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL"

EOM
    elif [[ $1 == "fitcons" ]];then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/wgsa_dbNSFP/resources/fitConsv1.01/fitcons_v1.01.txt.gz"
columns = [3]
names=["fitcons"]
ops=["self"]

EOM
    elif [[ $1 == "linsight" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/wgsa_dbNSFP/resources/LINSIGHT/LINSIGHT_hg19.bed.bgz"
columns=[4]
names=["linsight_g"]
ops=["self"]

EOM
    elif [[ $1 == "traP" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/traP/traP_v2.txt.gz"
columns=[5]
names=["traP"]
ops=["self"]

EOM
    elif [[ $1 == "ReMM" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/reMM_genomiser/ReMM.v0.3.1.txt.gz"
columns=[3]
names=["ReMM"]
ops=["self"]

EOM
    fi
}

tools=(fitcons linsight ReMM traP spliceai)
if [[ -f "$PWD/anno.conf" ]]; then
    rm "$PWD/anno.conf"
fi

if [ ! -z $4 ]; then
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
else
    for elem in "${tools[@]}"
        do
        annotate $elem
        done
fi

config=$(readlink -f $PWD/anno.conf)
if [ -f $PWD/custom.lua ]; then
    custom_lua=True
    lua_file="$PWD/custom.lua"
else
    custom_lua=False
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
#SBATCH --image=mcfonsecalab/variantutils:0.4

cd /home/pedro.barbosa/scratch/vep
mkdir \$SLURM_JOB_ID && cd \$SLURM_JOB_ID

if [ $custom_lua == "True" ]; then
    srun shifter vcfanno -p 1 -lua $lua_file $config $invcf > $(basename $outfile)
else
    srun shifter vcfanno -p 1 $config $invcf > $(basename $outfile)
fi
mv $(basename $outfile) $outdir
cd ../ && rm -rf \$SLURM_JOB_ID
EOL

sbatch $PWD/vcfAnno.sbatch


