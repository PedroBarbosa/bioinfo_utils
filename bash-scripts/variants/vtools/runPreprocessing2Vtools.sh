#!/bin/bash

display_usage(){
    printf "1st argument must be the filtered VCF.
2nd argument must be the list of samples to process.
3rd argument must be a flag to perform null imputation on input VCF. Values: [true|false].
4th argument must be a flag to indicate wheather samples should be renamed (pos|neg) samples. Values: [true|false]
5th argument must be the list of Sarcomeric genes.
6th argument is otpional. If set, gnomad population is treated as control. Values: [true|false]\n"
}

if [[ -z "$1" || -z "$2" || -z "$3" || -z "$4" || -z "$5" ]]; then
    printf "Please set the required arguments.\n"
    display_usage
    exit 1
fi

vcf=$(readlink -f $1)
samples=$(readlink -f $2)
sarcomeric_genes=$(readlink -f $5)

if [ "$3" == "true" ]; then
    shifter --image=mcfonsecalab/variantutils:0.4 bcftools +setGT $vcf  -- --target-gt ./. --new-gt 0 | bgzip > ${vcf/.vcf.gz/_nullImputation.vcf.gz}
    vcfile=${vcf/.vcf.gz/_nullImputation.vcf.gz}
elif [ "$3" == "false" ]; then
    vcfile=$vcf
else
    printf "Please set valid value for 3rd argument.\n"
    display_usage
    exit 1
fi
shifter --image=mcfonsecalab/variantutils:0.4 bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" $vcfile > v.txt
shifter --image=mcfonsecalab/variantutils:0.4 python /home/pedro.barbosa/git_repos/bioinfo_utils/python-scripts/vcf-tools/printVCFfields.py $vcfile SYMBOL | awk 'NR>2' > g.txt
paste -d "\t" v.txt g.txt | awk -v OFS="\t" '{print $0}' > gene.mapping.txt
cat gene.mapping.txt | awk -v OFS="\t" '{if ($5=="ACTC1" || $5=="MYBPC3" || $5=="MYH7" || $5=="MYL2" || $5=="MYL3" || $5=="TNNI3" || $5=="TNNT2" || $5=="TPM1" ){print $0,1} else {print $0,0}}' > gene.mapping2.txt
mv gene.mapping2.txt gene.mapping.txt
printf "Unique gene first consequences found:\n"
cut -f5 gene.mapping.txt | sort | uniq -c 
sed -i 's/SPI1/MYBPC3/g' gene.mapping.txt
rm v.txt g.txt

if [[ "$4" == "true" && -z "$6" ]];then
    shifter --image=mcfonsecalab/variantutils:0.4 python /home/pedro.barbosa/git_repos/bioinfo_utils/python-scripts/vcf-tools/inspectSampleGenotypes.py -n $vcfile genotypes
    fgrep -vf genotypes_samplesWithNoVariant.txt $samples | sed 's/^/ctrl_/'g > pos_withTag.txt
    mv genotypes_samplesWithNoVariant.txt neg_withTag.txt
    sed -i 's/^/case_/'g neg_withTag.txt
    sed -i -e '$a\' neg_withTag.txt
    cat neg_withTag.txt pos_withTag.txt | sort -t _ -k2 > samples.rename.txt
    shifter --image=mcfonsecalab/variantutils:0.4 bcftools reheader -s samples.rename.txt $vcfile > ${vcfile/.vcf.gz/_processed.vcf.gz}
    rm genotypes*
    cat neg_withTag.txt pos_withTag.txt | awk -v OFS="\t" '{if ($0 ~ /case/) {print $0,1} else {print $0,0}}' > phenotype.txt

elif [[ "$4" == "false" && -z "$6" || "$6" == "false" ]];then
    shifter --image=mcfonsecalab/variantutils:0.4 bcftools query -l $vcfile | awk -v OFS="\t" '{if ($0 ~ /case/) {print $0,1} else {print $0,0}}' > phenotype.txt
elif [[ "$6" == "true" ]]; then
    shifter --image=mcfonsecalab/variantutils:0.4 bcftools query -l $vcfile | awk -v OFS="\t" '{if ($0 ~ /ind/) {print $0,0} else {print $0,1}}' > phenotype.txt
else
    printf "Please set valid value for the 4th argument.\n"
    display_usage
    exit1
fi
sed -i '1 i\sample_name\tisPositive' phenotype.txt



