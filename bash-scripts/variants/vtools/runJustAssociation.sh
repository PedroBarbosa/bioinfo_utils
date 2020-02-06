#!/bin/bash

display_usage(){
    printf "
1th argument is the column from the phenotype file to test (e.g. isPositive)
2th argument is the column in the variants/gene mapping file to group the statistical tests (e.g. gene_symbol)
3th argument is the outdir
4th argument is the basename of the output files.
5th argument is optional. If set, it will discard samples and/or variants that are NA in more than the given fraction. (e.g. 0.1) \n"
}

if [[ -z "$1" || -z "$2" || -z "$3" || -z "$4" ]] ; then
    printf "Please set the required arguments for the script\n"
    display_usage
    exit 1
else
    phenotype=$1
    group=$2
    if [ ! -d $(readlink -f $3) ]; then 
        mkdir $(readlink -f $3)
    fi
    outdir=$(readlink -f $3)
    outbasename=$4
fi

CMD="vtools associate variant $phenotype --group_by $group --force -j8"
if [ ! -z "$5" ];then
    CMD="$CMD --discard_samples '%(NA)>$5' --discard_variants '%(NA)>$5'"
fi 

$CMD --method "CFisher --name CMC_fisher" > $outdir/${outbasename}_CMC_fisher.txt 
$CMD --method "BurdenBt --name Burden"    > $outdir/${outbasename}_burden.txt 
#$CMD --method "KBAC --name kbac -p 5000 --mafupper 0.2"  > $outdir/${outbasename}_KBAC.txt
#$CMD --method "Calpha --name Calpha --mafupper 0.2"      > $outdir/${outbasename}_Calpha.txt
#$CMD --method "RBT --name RBT -p 5000 --mafupper 0.2"    > $outdir/${outbasename}_rbt.txt 
#$CMD --method "VTtest --name vt -p 5000 --mafupper 0.2"  > $outdir/${outbasename}_vt.txt 
#$CMD --method "aSum --name aSum -p 5000 --mafupper 0.2"  > $outdir/${outbasename}_asum.txt

for file in $(find $outdir -name "${outbasename}*"); do awk 'NR == 1; NR > 1 {print $0 | "sort -V"}' $file > ${file/.txt/_sort.txt}; done

#paste -d "\t" $outdir/*sort.txt | cut -f1,2,3,4,6,14,21,27,33,41,49 > $outdir/final_${outbasename}.txt
paste -d "\t" $outdir/*sort.txt | cut -f1,2,3,4,6,13 > $outdir/final_${outbasename}.txt
#cat $outdir/final_${outbasename}.txt | vtools_report plot_association qq --label_top 8 -s --color Dark2 -b -o $outdir/plot_${outbasename}
cat ${file/.txt/_sort.txt} | vtools_report plot_association qq --label_top 8 -s --color Dark2 -b -o $outdir/plot_${outbasename}
