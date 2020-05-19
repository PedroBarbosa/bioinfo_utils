#!/bin/bash
display_usage(){
    printf "
1th argument is VCF file.
2rd argument is the list of genes to analyze.
3th argument is the file with the samples in the vcf.
4th argument refers to the flag whether gnomad is the control population. [true|false]. Default: false.
5th argument is the outdir.\n"
}

if [[ -z "$1" || -z "$2" || -z "$3" || -z "$4" || -z "$5" ]] ; then
    printf "Please set the required arguments for the script\n"
    display_usage
    exit 1
else
    vcf=$(readlink -f $1)
    listgenes=$(readlink -f $2)
    samples=$(readlink -f $3)
    if [ ! -d $(readlink -f $5) ]; then
        mkdir $(readlink -f $5)
    fi
    outdir=$(readlink -f $5)

    if [[ $4 == "false" ]];then
        cat $samples | awk -v OFS="\t" '{if ($0 ~ /case/) {print $0}}'  > $outdir/pos.txt
	cat $samples | awk -v OFS="\t" '{if ($0 ~ /control/) {print $0}}' > $outdir/neg.txt
    elif [[ $4 == "true" ]]; then
        cat $samples | awk -v OFS="\t" '{if ($0 ~ /H/) {print $0}}'  > $outdir/pos.txt
	cat $samples | awk -v OFS="\t" '{if ($0 ~ /ind/) {print $0}}'  > $outdir/neg.txt
    else
        printf "Please set a valid value for the 4th argument.\n"
        display_usage
        exir 1
    fi
    outdir=$(readlink -f $5)
fi
npos=$(cat $outdir/pos.txt | wc -l)
nneg=$(cat $outdir/neg.txt | wc -l)

vep_image="shifter --image=ensemblorg/ensembl-vep:latest"
bcftools_image="shifter --image=mcfonsecalab/variantutils:0.4 bcftools"
inspectSampleGenotypes="shifter --image=mcfonsecalab/variantutils:0.4 python /home/pedro.barbosa/git_repos/bioinfo_utils/python-scripts/vcf-tools/inspectSampleGenotypes.py"
tmp_remove="ls -la /tmp/ | grep 'pedro.barbosa' | awk ' { print \$9 } '"

cat > $PWD/variantCounts.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=variantCounts
#SBATCH --time=72:00:00
#SBATCH --mem=150G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --workdir=/home/pedro.barbosa/scratch/vep
#SBATCH --output=/home/pedro.barbosa/scratch/vep/variantCounts_%j.log
#SBATCH --image=ensemblorg/ensembl-vep:latest

mkdir \$SLURM_JOB_ID && cd \$SLURM_JOB_ID

while read line; do
    mkdir \$line
    srun shifter filter_vep -i $vcf --vcf_info_field ANN --only_matched -f "SYMBOL is \$line" | $bcftools_image view -c1 - > \$line/\$line.vcf
    echo "Gene: \$line"
    echo -e "Pos:$npos\nNeg:$nneg"
    echo "#Variants: \$($bcftools_image view -H \$line/\$line.vcf | wc -l)"
    $bcftools_image view -S $outdir/neg.txt \$line/\$line.vcf > \$line/neg.vcf
    $bcftools_image view -S $outdir/pos.txt \$line/\$line.vcf > \$line/pos.vcf
    $inspectSampleGenotypes -n \$line/pos.vcf \$line/pos
    $inspectSampleGenotypes -n \$line/neg.vcf \$line/neg
    npos_novariant=\$(cat \$line/pos_samplesWithNoVariant.txt | wc -l)
    nneg_novariant=\$(cat \$line/neg_samplesWithNoVariant.txt | wc -l)
    pos_withAnyVariant=\$((npos - npos_novariant))
    neg_withAnyVariant=\$((nneg - nneg_novariant))
    percent_pos=\$(awk "BEGIN { pc=\${pos_withAnyVariant}/${npos}*100; print pc}")
    percent_neg=\$(awk "BEGIN { pc=\${neg_withAnyVariant}/${nneg}*100; print pc}")
    echo "Num positives with any variant: \$pos_withAnyVariant (\$percent_pos%)"
    echo -e "Num negatives with any variant: \$neg_withAnyVariant (\$percent_neg%)\n"

done < $listgenes

mv * ../*\${SLURM_JOB_ID}.log $outdir
cd ../ && rm -rf \$SLURM_JOB_ID
rm -rf $tmp_remove
EOL
sbatch variantCounts.sbatch

