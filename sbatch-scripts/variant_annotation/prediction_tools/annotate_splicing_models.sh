#!/bin/bash
display_usage(){
    printf "
1st argument is the VCF file to annotate.
2nd argument is the output file. 
3rd argument is the output directory.
4th argument is the models to annotate, comma separated. Default: no model annotation.
5th argument are the result files of each model. Must come in the same order as the previous argument.
6th argument is optional. Do not reate and run a slurm sbatch script on the fly. Default: false, it creates. Values: [true|false|-]. Set '-' to skip the argument.\n

Possible tools:
HAL
mmsplice_deltaLogitPSI
mmsplice_pathogenicity
mmsplice_efficiency
kipoisplice_4
kipoisplice_4cons
maxentscan_5
maxentscan_3
SCAP
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

possible_models=(HAL mmsplice_deltaLogitPSI mmsplice_pathogenicity mmsplice_efficiency kipoisplice_4 kipoisplice_4cons maxentscan_5 maxentscan_3 SCAP)
annotate(){

    if [[ $1 == "HAL" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="$2"
fields=["KV:kipoi:HAL:DIFF"]
names=["HAL_DIFF"]
ops=["self"]

EOM

    elif [[ $1 == "mmsplice_deltaLogitPSI" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="$2"
columns=[5]
names=["mmsplice_deltaPSI"]
ops=["self"]

EOM
    
    elif [[ $1 == "mmsplice_pathogenicity" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]] 
file="$2"
columns=[5]
names=["mmsplice_patho"]
ops=["self"]

EOM

    elif [[ $1 == "mmsplice_efficiency" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="$2"
columns=[5]
names=["mmsplice_efficiency"]
ops=["self"]

EOM

    elif [[ $1 == "kipoisplice_4" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="$2"
columns=[5]
names=["kipoisplice_4"]
ops=["self"]

EOM


elif [[ $1 == "kipoisplice_4cons" ]]; then
    cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="$2"
columns=[5]
names=["kipoisplice_4cons"]
ops=["self"]

EOM

   elif [[ $1 == "maxentscan_5" ]]; then 
cat <<EOM >>$PWD/anno.conf 
[[annotation]]
file="$2"
fields = ["KV:kipoi:MaxEntScan/5prime:DIFF"]
names=["maxentscan_5"]
ops=["self"]

EOM

    elif [[ $1 == "maxentscan_3" ]]; then 
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="$2"
fields=["KV:kipoi:MaxEntScan/3prime:DIFF"]
names=["maxentscan_3"]
ops=["self"]

EOM

    elif [[ $1 == "SCAP" ]]; then 
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/S-CAP/scap_COMBINED_v1.0.vcf.gz"
fields=["SCAP","SCAP"]
names=["SCAP","SCAP"]
ops=["self","self"]

EOM
fi
}


if [[ -f "$PWD/anno.conf" ]]; then
    rm "$PWD/anno.conf"
fi


if [[ $4 == "-" || -z "$4" ]]; then
    printf "Don't you want to annotate your VCF with any model? Please set at least one.\n"
    display_usage
    exit 1
    
else
    IFS=','
    read -r -a array <<< "$4"
    read -r -a files <<< "$5"
    i=0
    for elem in "${array[@]}"
        do
            if [[ ! " ${possible_models[@]} " =~ " ${elem} " ]]; then
                printf "${elem} is not included in the list of available models"
                display_usage
                exit 1
            elif [[ "$elem" == "SCAP" ]]; then #if models with precomputed scores
                annotate $elem
            else
                fullpath_file=$(readlink -f ${files[$i]})
                annotate $elem $fullpath_file
                i=$((i + 1))
            fi
    done
fi

config=$(readlink -f $PWD/anno.conf)
cmd="vcfanno $config $invcf | shifter --image=ummidock/ubuntu_base:latest bgzip > $(basename $outfile)"

if [[ -z "$6" || $6 == "false" || $6 == "-" ]]; then
   cat > $PWD/vcfAnnoSplicing.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=vcfAnnotSplicing
#SBATCH --time=72:00:00
#SBATCH --mem=150G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
##SBATCH --workdir=/home/pedro.barbosa/scratch/vep
#SBATCH --output=%j_vcfannoSplicing.log
#SBATCH --image=mcfonsecalab/variantutils:0.5

cd /home/pedro.barbosa/scratch/vep
mkdir \${SLURM_JOB_ID}_vcfAnno && cd \${SLURM_JOB_ID}_vcfAnno
srun shifter $cmd
mv $(basename $outfile) $outdir
cd ../ && rm -rf \$SLURM_JOB_ID
EOL
    sbatch $PWD/vcfAnnoSplicing.sbatch

elif [[ "$6" == "true" ]];then
    echo "srun shifter --image=mcfonsecalab/variantutils:0.5 $cmd"
    shifter --image=mcfonsecalab/variantutils:0.5 $cmd
else
    printf "Please set a valid value for the 6th argument\n"
    display_usage
    exit 1
fi

