#!/bin/bash
display_usage(){
    printf "
1st argument is the input VCF.
2nd argument is the output basename. 
3rd argument is the output directory.
4th argument is optional. Refers to the models to run, comma separated. Default: all
5th argument is optional. Refers to the fasta file. Default: '/home/pedro.barbosa/mcfonseca/shared/genomes/human/hg19/GRCh37.primary_assembly_nochr.genome.fa'
6th argument is optional. Refers to the gtf file. Default: '/home/pedro.barbosa/mcfonseca/shared/genomes/human/hg19/gencode.v28lift37.annotation.nochr.gtf'

Possible models:
HAL
mmsplice_deltalogitPSI
mmsplice_pathogenicity
mmsplice_efficiency
kipoisplice_4
kipoisplice_4cons
maxentscan_5prime
maxentscan_3prime
"
}

if [[ -z $1 || -z $2 || -z $3 ]]; then
    printf "Error. Please set the required arguments for the script.\n"
    display_usage
    exit 1
fi
infile=$(readlink -f "$1")
if [[ ${infile: -4} == ".bgz" ]]; then
    zcat $infile > ${infile/.bgz/} 
    invcf=${infile/.bgz/}

elif [[ ${infile: -3} == ".gz" ]]; then
    zcat $infile > ${infile/.gz/}
    invcf=${infile/.gz/}
else
    invcf=$infile
fi
outbasename=$2
outdir=$(readlink -f $3)
outfinal="$outdir/$outbasename"
possible_models=(HAL kipoisplice_4 kipoisplice_4cons maxentscan_5prime maxentscan_3prime mmsplice_deltalogitPSI mmsplice_pathogenicity mmsplice_efficiency mmsplice)

base_sbatch (){
cat > $PWD/runKipoi.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=kipoi_predictions
#SBATCH --time=72:00:00
#SBATCH --mem=150G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=%j_kipoi_predictions.log
#SBATCH --image=mcfonsecalab/variantutils:0.5

cd /home/pedro.barbosa/scratch/vep
mkdir \${SLURM_JOB_ID}_kipoi && cd \${SLURM_JOB_ID}_kipoi
source activate kipoi-shared__envs__kipoi-py3-keras2

EOL

}

#$1 = vcf
#$2 = fasta
#$3 = gtf
#$4 = outbasename
HAL () {
cat >> $PWD/runKipoi.sbatch <<EOL
printf "`date` INFO: HAL started.\n"
kipoi veff score_variants HAL -o ${4}_HAL.vcf --dataloader_args='{"gtf_file":"$3", "fasta_file":"$2"}' -i "$1"
awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,2 -V "}' ${4}_HAL.vcf | bcftools norm -d none | shifter --image=ummidock/ubuntu_base:latest bgzip > ${4}_HAL.vcf.gz
shifter --image=ummidock/ubuntu_base:latest tabix -p vcf ${4}_HAL.vcf.gz
printf "`date` INFO: HAL finished.\n"

EOL
}

kipoiSplice_4 () {
cat >> $PWD/runKipoi.sbatch <<EOL
printf "`date` INFO: kipoiSplice_4 started.\n"
kipoi predict KipoiSplice/4 -o ${4}_kipoisplice_4.tsv --dataloader_args='{"fasta_file":"$2", "gtf_file":"$3", "vcf_file":"$1"}'
awk -v OFS="\t" '{ print \$2, \$4, \$5, \$1, \$7}' ${4}_kipoisplice_4.tsv | tail -n+2 | sort -k1,2 -V > ${4}_kipoisplice_4_to_annotate.tsv
sed  -i $'1i#chrom\tpos\tref\talt\tscore' ${4}_kipoisplice_4_to_annotate.tsv && shifter --image=ummidock/ubuntu_base:latest bgzip ${4}_kipoisplice_4_to_annotate.tsv
shifter --image=ummidock/ubuntu_base:latest tabix -s1 -b2 -e2 ${4}_kipoisplice_4_to_annotate.tsv.gz
printf "`date` INFO: kipoiSplice_4 finished.\n"

EOL
}

kipoiSplice_4cons () {
cat >> $PWD/runKipoi.sbatch <<EOL
printf "`date` INFO: kipoiSplice_4cons started.\n"
kipoi predict KipoiSplice/4cons -o ${4}_kipoisplice_4cons.tsv --dataloader_args='{"fasta_file":"$2", "gtf_file":"$3", "vcf_file":"$1"}'
awk -v OFS="\t" '{ print \$2, \$4, \$5, \$1, \$7}' ${4}_kipoisplice_4cons.tsv | tail -n+2 | sort -k1,2 -V > ${4}_kipoisplice_4cons_to_annotate.tsv
sed  -i $'1i#chrom\tpos\tref\talt\tscore' ${4}_kipoisplice_4cons_to_annotate.tsv && shifter --image=ummidock/ubuntu_base:latest bgzip ${4}_kipoisplice_4cons_to_annotate.tsv
shifter --image=ummidock/ubuntu_base:latest tabix -s1 -b2 -e2 ${4}_kipoisplice_4cons_to_annotate.tsv.gz
printf "`date` INFO: kipoiSplice_4cons finished.\n"

EOL
}

maxentscan_5 () {
cat >> $PWD/runKipoi.sbatch <<EOL
printf "`date` INFO: maxentscan_5 started.\n"
echo "aii ${4}"
kipoi veff score_variants MaxEntScan/5prime -o ${4}_maxentscan_5.vcf --dataloader_args='{"gtf_file":"$3", "fasta_file":"$2"}' -i "$1"
awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,2 -V "}' ${4}_maxentscan_5.vcf | bcftools norm -d none | shifter --image=ummidock/ubuntu_base:latest bgzip > ${4}_maxentscan_5.vcf.gz
shifter --image=ummidock/ubuntu_base:latest tabix -p vcf ${4}_maxentscan_5.vcf.gz
printf "`date` INFO: maxentscan_5 finished.\n"

EOL
}

maxentscan_3 () {
cat >> $PWD/runKipoi.sbatch <<EOL
printf "`date` INFO: maxentscan_3 started.\n"
kipoi veff score_variants MaxEntScan/3prime -o ${4}_maxentscan_3.vcf --dataloader_args='{"gtf_file":"$3", "fasta_file":"$2"}' -i "$1"
awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,2 -V "}' ${4}_maxentscan_3.vcf | bcftools norm -d none | shifter --image=ummidock/ubuntu_base:latest bgzip > ${4}_maxentscan_3.vcf.gz
shifter --image=ummidock/ubuntu_base:latest tabix -p vcf ${4}_maxentscan_3.vcf.gz
printf "`date` INFO: maxentscan_3 finished.\n"

EOL
}

mmsplice (){
cat >> $PWD/runKipoi.sbatch <<EOL
conda deactivate
source activate kipoi-MMSplice
printf "`date` INFO: mmsplice delta logit PSI started.\n"
kipoi predict MMSplice/deltaLogitPSI -o ${4}_mmsplice_deltaLogitPSI.tsv --dataloader_args='{"fasta_file":"$2", "gtf":"$3", "vcf_file":"$1"}'
awk -v OFS="\t" '{ print \$21, \$23, \$24, \$20, \$26}' ${4}_mmsplice_deltaLogitPSI.tsv | tail -n+2 > ${4}_mmsplice_deltaLogitPSI_to_annotate.tsv
sed -i $'1i#chrom\tpos\tref\talt\tscore' ${4}_mmsplice_deltaLogitPSI_to_annotate.tsv && shifter --image=ummidock/ubuntu_base:latest bgzip  ${4}_mmsplice_deltaLogitPSI_to_annotate.tsv
shifter --image=ummidock/ubuntu_base:latest tabix -s1 -b2 -e2 ${4}_mmsplice_deltaLogitPSI_to_annotate.tsv.gz
printf "`date` INFO: mmsplice delta logit PSI finished.\n"

printf "`date` INFO: mmsplice efficiency started.\n"
kipoi predict MMSplice/deltaLogitPSI -o ${4}_mmsplice_efficiency.tsv --dataloader_args='{"fasta_file":"$2", "gtf":"$3", "vcf_file":"$1"}'
awk -v OFS="\t" '{ print \$21, \$23, \$24, \$20, \$26}' ${4}_mmsplice_efficiency.tsv | tail -n+2 > ${4}_mmsplice_efficiency_to_annotate.tsv
sed -i $'1i#chrom\tpos\tref\talt\tscore' ${4}_mmsplice_efficiency_to_annotate.tsv && shifter --image=ummidock/ubuntu_base:latest bgzip ${4}_mmsplice_efficiency_to_annotate.tsv
shifter --image=ummidock/ubuntu_base:latest tabix -s1 -b2 -e2 ${4}_mmsplice_efficiency_to_annotate.tsv.gz
printf "`date` INFO: mmsplice efficiency finished.\n"

printf "`date` INFO: mmsplice pathogenicity started.\n"
kipoi predict MMSplice/pathogenicity -o ${4}_mmsplice_pathogenicity.tsv --dataloader_args='{"fasta_file":"$2", "gtf":"$3", "vcf_file":"$1"}'
awk -v OFS="\t" '{ print \$21, \$23, \$24, \$20, \$27}' ${4}_mmsplice_pathogenicity.tsv | tail -n+2 > ${4}_mmsplice_pathogenicity_to_annotate.tsv
sed -i $'1i#chrom\tpos\tref\talt\tscore' ${4}_mmsplice_pathogenicity_to_annotate.tsv && shifter --image=ummidock/ubuntu_base:latest bgzip ${4}_mmsplice_pathogenicity_to_annotate.tsv
shifter --image=ummidock/ubuntu_base:latest tabix -s1 -b2 -e2 ${4}_mmsplice_pathogenicity_to_annotate.tsv.gz
printf "`date` INFO: mmsplice pathogenicity finished.\n"
EOL
} 

base_sbatch

if [[ -z $5 || "$5" == "-" ]]; then
    fasta="/home/pedro.barbosa/mcfonseca/shared/genomes/human/hg19/GRCh37.primary_assembly_nochr.genome.fa"
else
    fasta=$(readlink -f "$5") 
fi

if [[ -z $6 || "$6" == "-" ]]; then
    gtf="/home/pedro.barbosa/mcfonseca/shared/genomes/human/hg19/gencode.v28lift37.annotation.nochr.gtf"
else
    gtf=$(readlink -f "$6")
fi

if [[ $4 == "all" || $4 == "-" || -z "$4" ]]; then
    printf "All models will be run by default.\n"
    HAL $invcf $fasta $gtf $outfinal
    kipoiSplice_4 $invcf $fasta $gtf $outfinal 
    maxentscan_5 $invcf $fasta $gtf $outfinal
    maxentscan_3 $invcf $fasta $gtf $outfinal
    mmsplice $invcf $fasta $gtf $outfinal
else
    IFS=','
    read -r -a array <<< "$4"
    for elem in "${array[@]}";
    do
            if [[ ! " ${possible_models[@]} " =~ " ${elem} " ]]; then
                printf "${elem} is not included in the list of available models"
                display_usage
                exit 1
            else
                if [[ "$elem" == "HAL" ]]; then
                    HAL $invcf $fasta $gtf $outfinal
                elif [[ "$elem" == "kipoisplice_4" ]]; then
                    kipoiSplice_4 $invcf $fasta $gtf $outfinal
                elif [[ "$elem" == "kipoisplice_4cons" ]];then
                    kipoiSplice_4cons $invcf $fasta $gtf $outfinal
                elif [[ "$elem" == "maxentscan_5prime" ]]; then
                    maxentscan_5 $invcf $fasta $gtf $outfinal
		elif [[ "$elem" == "maxentscan_3prime" ]]; then
                    maxentscan_3 $invcf $fasta $gtf $outfinal
                elif [[ "$elem" == "mmsplice" || "$elem" == "mmsplice_deltalogitPSI" || "$elem" == "mmsplice_pathogenicity" || "$elem" == "mmsplice_efficiency" ]]; then
                    mmsplice $invcf $fasta $gtf $outfinal
                fi
            fi
    done
fi

sbatch $PWD/runKipoi.sbatch

