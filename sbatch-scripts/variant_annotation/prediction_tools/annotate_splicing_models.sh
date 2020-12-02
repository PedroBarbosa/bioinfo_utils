#!/bin/bash
display_usage(){
    printf "
1st argument is the VCF file to annotate.
2nd argument is the output file. (output will be automatically gzipped) 
3rd argument is the output directory.
4h argument is the models to annotate, comma separated. Default: no model annotation.
5th argument are the result files of each model. Must come in the same order as the previous argument.\n

Possible tools:
HAL
mmsplice_deltaLogitPSI
mmsplice_pathogenicity
mmsplice_efficiency
kipoisplice_4
kipoisplice_4cons
maxentscan_5
maxentscan_3
rbp_clip/TARGET_NAME
"
}

if [[ -z $1 || -z $2 || -z $3 ]]; then
    printf "Error. Please set the required arguments for the script.\n"
    display_usage
    exit 1
fi
invcf=$(readlink -f "$1")
outfile=$(readlink -f "$2")
if [[ -z "$outfile" ]];then
    printf "Error on the output file. Please set a different string. Perhaps remove the gz extension at the end?\n"
    display_usage
    exit 1
fi
outdir=$(readlink -f $3)
possible_models=(HAL mmsplice_deltaLogitPSI mmsplice_pathogenicity mmsplice_efficiency kipoisplice_4 kipoisplice_4cons maxentscan_5 maxentscan_3 rbp_eclip)
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

    elif [[ $1 == "rbp_eclip" ]]; then
cat <<EOM >>$PWD/anno.conf
[[annotation]]
file="$2"
fields=["KV:kipoi:$3:DIFF"]
names=["$3"]
ops=["self"]

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
                if [[ "${elem}"  == *"rbp_eclip"* ]]; then
                    fullpath_file=$(readlink -f ${files[$i]})
                    annotate rbp_eclip $fullpath_file ${elem}
                    i=$((i + 1))
                else
                    printf "${elem} is not included in the list of available models"
                    display_usage
                    exit 1
                fi
            else
                
                fullpath_file=$(readlink -f ${files[$i]})
                annotate $elem $fullpath_file
                i=$((i + 1))
            fi
    done
fi

config=$(readlink -f $PWD/anno.conf)
cmd="vcfanno $config $invcf | shifter --image=ummidock/ubuntu_base:latest bgzip > ${outfile}"

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
#SBATCH --image=mcfonsecalab/variantutils:latest

echo "$cmd"
srun shifter $cmd
srun shifter --image=ummidock/ubuntu_base:latest tabix -p vcf ${outfile}
if [ ! -f "${outdir}/$(basename $outfile)"  ];then
    echo "${outdir}/$(basename $outfile)"
    mv ${outfile}* $outdir
fi

cd ../ && rm -rf \$SLURM_JOB_ID
EOL
sbatch $PWD/vcfAnnoSplicing.sbatch

