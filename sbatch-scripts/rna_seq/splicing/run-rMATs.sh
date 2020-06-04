#!/usr/bin/env bash
display_usage() {
echo 'Script to run rMATs for multiple BAM files.

-1st argument must be the file listing the replicates for condition 1. Replicates must be split by ",".
-2st argument must be the file listing the replicates for condition 2. Replicates must be split by ","
-3rd argument must be the GTF annotation file (must be unzipped). Use "-" to skip argument and use default. Default: gencode hg38 v33 primary assembly.
-4th argument must be the output directory. (will serve also as basename of main output files)
-5th argument must be the labels for each group, split by ",".
-6th argument is optional. Refers to the read type. Values [paired|single|-]. Default: paired.
-7th argument is optional. Refers to the library type of the RNA. Values: [fr-firstsrand|fr-secondstrand|fr-unstranded|-]. Default: fr-firststrand.
-8th argument is optional. Refers to whether statistical analysis should be performed. Values: [true|false|-]. Default: true
-9th argument is optional. Refers to the read length. Default: not set.
-10th argument is optional. Refers to the threshold for significance to apply, each separated by comma (min_avg_reads, fdr_threshold, deltaPSI_threshold). Default:5,0.05,0.2
-11th argument is optional. If a previous rmats run was performed, pipeline will start from output files stored in the output directory (4th argument). Values:[true|false].Default:false'

}


if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] ; then
        printf "ERROR:Please provide at least the 4 first arguments required for the script.\n\n"
        display_usage 
        exit 1
fi

BAM_1=$(readlink -f "$1")
BAM_2=$(readlink -f "$2")
if [[ "$3" == "-" ]]; then
    GTF="/home/pedro.barbosa/mcfonseca/shared/genomes/human/hg38/gencode.v33.primary_assembly.annotation.gtf"
else
    GTF=$(readlink -f "$3")
fi
if [[ ! -d $(readlink -f "$4") ]]; then
    mkdir $(readlink -f "$4")
fi
OUTDIR=$(readlink -f "$4")
OUT_BASENAME=$4

IFS=','
read -r -a array <<< "$5"
label1=${array[0]} 
label2=${array[1]}

if [[ -z "$6" || "$6" == "paired" || "$6" == "-" ]]; then
    READTYPE="paired"
elif [[ "$6" == "single" ]]; then
    READTYPE="single"
else
    printf "Please set a valid value for 6th arg.\n"
    display_usage
    exit 1
fi

if [[ -z "$7" || "$7" == "fr-firststrand" || "$7" == "-" ]]; then
    LIBTYPE="fr-firststrand"
elif [[ "$7" == "fr-secondstrand" ]]; then
    LIBTYPE="$7"
elif [[ "$7" == "fr-unstranded" ]]; then
    LIBTYPE="$7"
else
    printf "Please set a valid value for 7th arg.\n"
    display_usage
    exit 1
fi

if [[ -z "$8" || "$8" == "true" || "$8" == "-" ]]; then
    DO_STAT_ANALYSIS="true"
    STATS="true"
elif [[ "$8" == "false" ]]; then
    DO_STAT_ANALYSIS="false"
    STATS="false"
else
    printf "Please set a valid value for 8th arg.\n"
    display_usage
    exit 1
fi 

if [[ -z "$9" || "$9" == "-" ]]; then
    CMD="rmats4.py --b1 $BAM_1 --b2 $BAM_2 --gtf $GTF --od \$PWD -t paired --nthread \$SLURM_CPUS_PER_TASK --cstat 0.0001 --anchorLength 3"
    #CMD="Darts_BHT rmats_count --b1 $BAM_1 --b2 $BAM_2 --gtf $GTF --od \$PWD -t paired --nthread \$SLURM_CPUS_PER_TASK"
else
    CMD="rmats4.py --b1 $BAM_1 --b2 $BAM_2 --gtf $GTF --od \$PWD -t paired --nthread \$SLURM_CPUS_PER_TASK --cstat 0.0001 --readLength $9 --anchorLength 3"
fi
#--enable-unicode=ucs4

if [[ $LIBTYPE == "paired" ]];then
    CMD="$CMD -libType $LIBTYPE"
fi

if [[ $DO_STAT_ANALYSIS == "false" ]]; then
    CMD="$CMD --statoff"
fi


if [[ -z "${10}" || "${10}" == "-" ]]; then
    min_avg_reads=5
    fdr_threshold=0.05
    deltaPSI_threshold=0.2
else
    IFS=","
    read -r -a array <<< "${10}"
    min_avg_reads=${array[0]}
    fdr_threshold=${array[1]}
    deltaPSI_threshold=${array[2]}
fi

if [[ -z "${11}" || "${11}" == "false" ]]; then
    previous_run="false"
elif [[ "${11}" == "true" ]]; then
    previous_run="true"
else
    printf "Please set a valid value for the 10th argument.\n"
    display_usage
    exit 1
fi

cat > rmats.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=rmats
#SBATCH --time=72:00:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=mcfonsecalab/rmats:latest
#SBATCH --output=%j_rmats.log

if [[ $previous_run == "false" ]];then
    workdir="/home/pedro.barbosa/scratch/rna_seq/splicing/\$SLURM_JOB_ID"
    mkdir \$workdir && cd \$workdir
    srun shifter $CMD
elif [[ $previous_run == "true" ]];then
    cd $OUTDIR
fi

if [[ $STATS == "true" ]]; then
    printf "Producing tables of significant events\n"
    printf "Tresholds used for the filtering (fdr, dPSI, min_avg_reads): $fdr_threshold, $deltaPSI_threshold, $min_avg_reads\n"
    srun shifter Rscript /python_env/run_maser.R $label1 $label2 $min_avg_reads $fdr_threshold $deltaPSI_threshold "JC" $GTF
    srun cut -f3 coverage_filt/* | sort | uniq | grep -v "geneSymbol" | sed 's/\\(ENSG[0-9]*\\)\\.[0-9]*/\\1/g' > ${OUT_BASENAME}_rmats_negative_list_to_GO.txt
    srun cut -f2 sign_events_* | sed 's/\\(ENSG[0-9]*\\)\\.[0-9]*/\\1/g' | sort | uniq > ${OUT_BASENAME}_rmats_positive_list_to_GO.txt
    srun echo "gene_id\tFDR\tdPSI\n" > ${OUT_BASENAME}_rmats_gene_ranks_to_GSEA.txt
    srun awk '{print \$3, \$(NF-3), \$NF}' OFS="\t" coverage_filt/* |  sed 's/\\(ENSG[0-9]*\\)\\.[0-9]*/\\1/g' | grep -v geneSymbol >> ${OUT_BASENAME}_rmats_gene_ranks_to_GSEA.txt

    printf "Done\nGenerating sashimi plots from all significant events..\n"
    events=("SE" "A5SS" "A3SS" "RI" "MXE")
    for e in "\${events[@]}"; do
        printf "Looking at significant \${e} events..\n"
        awk -F'\t' 'NR==FNR{c[\$1]++;next};c[\$1]' sign_events_\${e}.tsv \${e}.MATS.JC.txt > \${e}_toSashimi.txt
        mkdir \${e}_sashimi_plots
        CMD_SASHIMI="rmats2sashimiplot --b1 \$(cat $BAM_1) --b2 \$(cat $BAM_2) -t \${e} -e \${e}_toSashimi.txt --l1 $label1 --l2 $label2 --exon_s 1 --intron_s 5 -o \${e}_sashimi_plots"
        srun shifter \$CMD_SASHIMI
        cd \${e}_sashimi_plots
        mv Sashimi_index/* .
        mv Sashimi_plot/* .
        rm -rf *index* && rm -rf Sashimi_plot  
        cd ../ 
    done
else
    printf "Done. As stats were set to false. Most of the analysis won't be performed.\n"
fi

rename 's/sign/${OUT_BASENAME}_sign/g' *
mv * $OUTDIR
cd ../ && rm -rf \$SLURM_JOB_ID
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
EOL

sbatch rmats.sbatch
