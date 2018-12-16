#!/bin/bash
display_usage(){
    printf "
1th argument is the full directory where VCF is located (to mount within the VEP image with shifter)
2th argument is VCF file-
3rd argument is the list of genes.
4th argument is the file with the samples renamed. (pos and neg)
4th argument is the outdir.\n"
}

if [[ -z "$1" || -z "$2" || -z "$3"  ]] ; then
    printf "Please set the required arguments for the script\n"
    display_usage
    exit 1
else
    mountdir=$1
    vcf=$(basename $2)
    listgenes=$(readlink -f $3)
    samples_renamed=$(readlink -f $4)
    if [ ! -d $(readlink -f $5) ]; then 
        mkdir $(readlink -f $5)
    fi
    outdir=$(readlink -f $5)
fi
cat $samples_renamed | awk -v OFS="\t" '{if ($0 ~ /case/) {print $0}}'  > $outdir/pos.txt
cat $samples_renamed | awk -v OFS="\t" '{if ($0 ~ /control/) {print $0}}' > $outdir/neg.txt
npos=$(cat $outdir/pos.txt | wc -l)
nneg=$(cat $outdir/neg.txt | wc -l)

vep_image="shifter --image=ensemblorg/ensembl-vep:latest"
bcftools_image="shifter --image=mcfonsecalab/variantutils:0.4 bcftools"
inspectSampleGenotypes="shifter --image=mcfonsecalab/variantutils:0.4 python /home/pedro.barbosa/git_repos/bioinfo_utils/python-scripts/vcf-tools/inspectSampleGenotypes.py"
while read line; do
    mkdir $outdir/$line
    cd $outdir
    $vep_image -V $mountdir:/media filter_vep -i /media/$vcf --vcf_info_field ANN --only_matched -f "SYMBOL is $line" | $bcftools_image view -c1 - > $line/$line.vcf
    echo "Gene: $line"
    echo -e "Pos:$npos\nNeg:$nneg"
    echo "#Variants: $($bcftools_image view -H $line/$line.vcf | wc -l)"
    $bcftools_image view -S neg.txt $line/$line.vcf > $line/neg.vcf
    $bcftools_image view -S pos.txt $line/$line.vcf > $line/pos.vcf
    $inspectSampleGenotypes -n $line/pos.vcf $line/pos
    $inspectSampleGenotypes -n $line/neg.vcf $line/neg
    npos_novariant=$(cat $line/pos_samplesWithNoVariant.txt | wc -l)
    nneg_novariant=$(cat $line/neg_samplesWithNoVariant.txt | wc -l)
    pos_withAnyVariant=$((npos - npos_novariant))
    neg_withAnyVariant=$((nneg - nneg_novariant))
    percent_pos=$(awk "BEGIN { pc=${pos_withAnyVariant}/${npos}*100; print pc}")
    percent_neg=$(awk "BEGIN { pc=${neg_withAnyVariant}/${nneg}*100; print pc}")
    echo "Num positives with any variant: $pos_withAnyVariant ($percent_pos%)"
    echo -e "Num negatives with any variant: $neg_withAnyVariant ($percent_neg%)\n"

done < $listgenes
