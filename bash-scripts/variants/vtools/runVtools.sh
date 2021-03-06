#!/bin/bash

display_usage(){
    printf "1st argument is the project name.
2nd argument is the vcf file to import.
3rd argument is the phenotype file to import.
4th argument is the annotation database to import (*ann file).
5th argument is the variants/gene mapping file to import.
6th argument is the column from the phenotype file to test (e.g. isPositive)
7th argument is the column in the variants/gene mapping file to group the statistical tests (e.g. gene_symbol)
8th argument is the outdir
9th argument is the basename of the output files.
10th argument is the basename of the final output files.
11th argument is optional. If set, it will discard samples and/or variants that are NA in more than the given fraction. (e.g. 0.1)
\n"
}

if [[ -z "$1" || -z "$2" || -z "$3" || -z "$4" || -z "$5" || -z "$6" || -z "$7" || -z "$8" || -z "$9" || -z "${10}" ]] ; then
    printf "Please set the required arguments for the script\n"
    display_usage
    exit 1
else
    project=$1
    vcf=$(readlink -f $2)
    phenotype_file=$(readlink -f $3)
    ann_database=$(readlink -f $4)
    gene_mapping=$(readlink -f $5)
    phenotype=$6
    group=$7
   if [ ! -d $(readlink -f $8) ]; then 
        mkdir $(readlink -f $8)
    fi
    outdir=$(readlink -f $8)
    outbasename=$9
    finaloutbasename=${10}
fi

vtools remove project
vtools init $project
vtools import $vcf --build hg19 --var_info DP AC AN --geno_info DP_geno GQ_geno AD_geno
vtools phenotype --from_file $phenotype_file $phenotype
vtools use $ann_database -f $gene_mapping

#SINGLE VARIANT TESTS"
vtools update variant --from_stat 'num_gt_case=#(GT)' 'num_var_alleles_case=#(alt)' --samples "sample_name like 'H%'"
vtools update variant --from_stat 'num_gt_ctrl=#(GT)' 'num_var_alleles_ctrl=#(alt)' --samples "sample_name like 'ind%'"
vtools update variant --set "prop_pval=Fisher_exact(num_var_alleles_case, num_var_alleles_ctrl, 2*num_gt_case, 2*num_gt_ctrl)"
vtools output variant chr pos ref alt prop_pval | sort -r -k5 > $outdir/single_variant_fisher.txt

CMD="vtools associate variant $phenotype --group_by $group --force -j8"
if [ ! -z "${11}" ];then
    CMD="$CMD --discard_samples '%(NA)>${11}' --discard_variants '%(NA)>${11}'"
fi
$CMD --method "CFisher --name CMC_fisher --mafupper 0.2" > $outdir/${outbasename}_CMC_fisher.txt
#$CMD --method "BurdenBt --name Burden"    > $outdir/${outbasename}_burden.txt
$CMD --method "KBAC --name kbac -p 5000 --mafupper 0.2"  > $outdir/${outbasename}_KBAC.txt
#$CMD --method "Calpha --name Calpha --mafupper 0.2"      > $outdir/${outbasename}_Calpha.txt
#$CMD --method "RBT --name RBT -p 5000 --mafupper 0.2"    > $outdir/${outbasename}_rbt.txt
#$CMD --method "VTtest --name vt -p 5000 --mafupper 0.2"  > $outdir/${outbasename}_vt.txt
#$CMD --method "aSum --name aSum -p 5000 --mafupper 0.2"  > $outdir/${outbasename}_asum.txt

for file in $(find $outdir -name "${outbasename}*"); do awk 'NR == 1; NR > 1 {print $0 | "sort -V"}' $file > ${file/.txt/_sort.txt}; done

#paste -d "\t" $outdir/*sort.txt | cut -f1,2,3,4,6,14,21,27,33,41,49 > $outdir/final_${outbasename}.txt
paste -d "\t" $outdir/*sort.txt | cut -f1,2,3,4,6,12 > $outdir/final_${outbasename}.txt
if [ $group == "gene_symbol" ]; then
    nlabels=8
elif [ $group == "isSarcomeric" ];then
    nlabels=2
fi

cat $outdir/final_${outbasename}.txt | vtools_report plot_association qq --label_top $nlabels -s --color Dark2 -b -o $outdir/plot_${outbasename}
b1=$(echo ${finaloutbasename} | cut -f1 -d "_")
b2=$(echo ${finaloutbasename} | cut -f2 -d "_" | sed 's/\.\///g')
mv $outdir/final_${outbasename}.txt $outdir/${b1}_pvalues_${b2}.txt
mv $outdir/plot_${outbasename}.pdf $outdir/${b1}_plot_${b2}.pdf
cp $outdir/${b1}* ../

