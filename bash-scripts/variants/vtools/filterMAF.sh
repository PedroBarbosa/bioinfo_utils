#!/bin/bash

display_usage(){
    printf "
1th argument is the path where VCF is to mount image (pwd -P)
2nd argument is the VCF file
3th argument is the maf value.
4th argument is optional. If set, refers to the population (vep vcf field).
5th argument is optional. If set, it refers to the name of the output file.\n"
}
if [[ -z "$1" || -z "$2"  || -z "$3" ]] ; then
    printf "Please set the required arguments for the script\n"
    display_usage
    exit 1
else
    path=$1
    vcf=$(basename $2)
    if [[ $3 == *"./"* ]];then
        maf=$(echo "$3" | cut -f2 -d "/")
    else
        maf=$3
    fi

    if [[ -z "$4" ]]; then
        filter="(gnomADg_AF < $maf or not gnomADg_AF) and (gnomADg_AF_nfe < $maf or not gnomADg_AF_nfe) and (MAX_AF < $maf or not MAX_AF)"
    else
        filter="($4 < $maf or not $4) and (MAX_AF < $maf or not MAX_AF)"
    fi

    if [[ -z "$5" && ! -d "$maf" ]]; then
        mkdir $maf
        sleep 2
    elif [[ -z "$5" ]];then
        outfile="$maf/maffilt.vcf.gz"

    else
        outfile=$5
    fi
fi
vep_image="shifter --image=ensemblorg/ensembl-vep:latest"
bgzip="shifter --image=ummidock/ubuntu_base:latest bgzip"
$vep_image -V $path:/media filter_vep -i /media/$vcf --vcf_info_field ANN -f "$filter" | $bgzip > $outfile
