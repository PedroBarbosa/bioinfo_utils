#!/bin/bash

display_usage(){
    printf "
1st argument is the VCF file
2nd argument is the maf value.
3rd argument is optional. If set, refers to the population (vep vcf field). '-' skips its use.
4th argument is optional. If set, it refers to the name of the output file.\n"
}
if [[ -z "$1" || -z "$2"  ]] ; then
    printf "Please set the required arguments for the script\n"
    display_usage
    exit 1
else
    vcf=$(readlink -f $1)
    if [[ $2 == *"./"* ]];then
        maf=$(echo "$2" | cut -f2 -d "/")
    else
        maf=$2
    fi
    if [[ -z "$3" || "$3" == "-" ]]; then
        filter="(gnomADg_AF < $maf or not gnomADg_AF) and (gnomADg_AF_nfe < $maf or not gnomADg_AF_nfe) and (MAX_AF < $maf or not MAX_AF)"
    else
        filter="($3 < $maf or not $3) and (MAX_AF < $maf or not MAX_AF)"
    fi

    if [[ -z "$4" && ! -d "$maf" ]]; then
        mkdir $maf
        sleep 2
    elif [[ -z "$4" ]];then
        outfile="maffilt.vcf.gz"
    else
        outfile=$4
    fi
fi

if [[ -d "$maf" ]]; then
    cd $maf && outdir=$PWD && cd ../
else
    outdir=$PWD
fi

printf "Filter: $filter\n"
bgzip="shifter --image=ummidock/ubuntu_base:latest bgzip"
tmp_remove="ls -la /tmp/ | grep 'pedro.barbosa' | awk ' { print \$9 } '"
cat > $PWD/filter_byMAF_VEP.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=filterMAF_VEP
#SBATCH --time=72:00:00
#SBATCH --mem=150G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --workdir=/home/pedro.barbosa/scratch/vep
#SBATCH --output=/home/pedro.barbosa/scratch/vep/filterMAF_VEP_%j.log
#SBATCH --image=ensemblorg/ensembl-vep:latest

mkdir \$SLURM_JOB_ID && cd \$SLURM_JOB_ID
srun shifter filter_vep -i $vcf --only_matched --vcf_info_field ANN -f "$filter" | $bgzip > $outfile
mv $outfile ../*\${SLURM_JOB_ID}.log $outdir
cd ../ && rm -rf \$SLURM_JOB_ID
rm -rf $tmp_remove
EOL

sbatch filter_byMAF_VEP.sbatch
