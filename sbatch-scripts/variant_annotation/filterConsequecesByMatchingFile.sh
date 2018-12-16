#!/bin/bash
display_usage(){
    printf "
1st argument is the VCF file
2nd argument is VEP ANN field name to filter.
3rd argument is the file with the terms to filter.
4th argument is the name of the ouptut file.\n"
}

if [[ -z "$1" || -z "$2"  || -z "$3" || -z "$4" ]] ; then
    printf "Please set the required arguments for the script\n"
    display_usage
    exit 1
else
    vcf=$(readlink -f $1)
    ann_field=$2
    filter_file=$(readlink -f $3)
    out=$(readlink -f $4)
    outdir=$(dirname $out)
    outfile=$(basename $out)
    filter="$ann_field in $filter_file"
fi
bgzip="shifter --image=ummidock/ubuntu_base:latest bgzip"
tmp_remove="ls -la /tmp/ | grep 'pedro.barbosa' | awk ' { print \$9 } '"
cat > $PWD/filterVEP_byAuxFile.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=filterVEP_byAuxFile
#SBATCH --time=72:00:00
#SBATCH --mem=150G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --workdir=/home/pedro.barbosa/scratch/vep
#SBATCH --output=/home/pedro.barbosa/scratch/vep/filterVEP_fromAuxFile_%j.log
#SBATCH --image=ensemblorg/ensembl-vep:latest

mkdir \$SLURM_JOB_ID && cd \$SLURM_JOB_ID
srun shifter filter_vep -i $vcf --only_matched --vcf_info_field ANN -f "$filter" | $bgzip > $outfile
mv $outfile ../*\${SLURM_JOB_ID}.log $outdir
cd ../ && rm -rf \$SLURM_JOB_ID
rm -rf $tmp_remove
EOL

sbatch filterVEP_byAuxFile.sbatch
