#!/usr/bin/env bash
display_usage() {
echo 'Script to run rMATs for multiple BAM files.

-1st argument must be the file listing the replicates for condition 1. Replicates must be split by ",".
-2st argument must be the file listing the replicates for condition 2. Replicates must be split by ",". If not set (by using "-") PSI quantifications will be performed for the samples present in the 1st arg. No differential splicing will be tested. This is useful to have inclusion tables and perform PCA on the data, for instance.
-3rd argument must be the GTF annotation file (must be unzipped). Use "-" to skip argument and use default. Default: gencode hg38 v34 primary assembly.
-4th argument must be the output directory. (will serve also as basename of main output files).
-5th argument must be the labels for each group, split by ",". If 2nd argument is set to "-", this argument is irrelevant as no statistical comparisons between groups is done.
-6th argument is optional. Refers to whether statistical analysis should be performed. Values: [true|false|-]. Default: true
-7th argument is optional. Whether samples are paired and this should be considered in the statistical analysis. Values: [true|false|-]. Default: false
-8th argument is optional. If a previous rmats run was performed (--task prep or --task both), pipeline will start from rmats file (if set to "post") or directly to maser (if set to "maser") stored in the output directory (4th argument). Values:[false|post|maser|-].Default:false. If "-" is set, argument is skipped and defaults apply.
-9th argument is optional. Refers to the threshold for significance to apply, each separated by comma (min_avg_reads, fdr_threshold, deltaPSI_threshold). Default:5,0.05,0.2
-10th argument is optional. Refers to the read type. Values [paired|single|-]. Default: paired.
-11th argument is optional. Refers to the library type of the RNA. Values: [fr-firstsrand|fr-secondstrand|fr-unstranded|-]. Default: fr-firststrand.
-12th argument is optional. Refers to the read length. Default: 150bp.
'
}


if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] ; then
        printf "ERROR:Please provide at least the 4 first arguments required for the script.\n\n"
        display_usage 
        exit 1
fi

BAM_1=$(readlink -f "$1")
if [[ "$2" != "-" ]]; then
    ONLY_PSI="false"
    BAM_2=$(readlink -f "$2")
else
    ONLY_PSI="true"
fi

if [[ "$3" == "-" ]]; then
    GTF="/home/pedro.barbosa/mcfonseca/shared/genomes/human/hg38/gencode_v34_primary_assembly.gtf"
else
    GTF=$(readlink -f "$3")
fi
if [[ ! -d $(readlink -f "$4") ]]; then
    mkdir $(readlink -f "$4")
fi
OUTDIR=$(readlink -f "$4")
OUT_BASENAME=$4
CMD="rmats.py --gtf $GTF --od \$PWD -t paired --nthread \$SLURM_CPUS_PER_TASK --cstat 0.0001 --anchorLength 3 --novelSS --b1 $BAM_1"
#CMD="Darts_BHT rmats_count --b1 $BAM_1 --b2 $BAM_2 --gtf $GTF --od \$PWD -t paired --nthread \$SLURM_CPUS_PER_TASK"

# If two groups should be compared
if [[ $ONLY_PSI == "false" ]]; then
    CMD+=" --b2 $BAM_2"
    IFS=','
    read -r -a array <<< "$5"
    label1=${array[0]} 
    label2=${array[1]}
fi

#Statistical analysis
if [[ -z "$6" || "$6" == "true" || "$6" == "-" ]]; then
    if [[ "$ONLY_PSI" == "true" ]]; then
        printf "Statistical analysis (6th arg) should be disabled when aiming to quantify PSIs (2nd argument set to '-').\n"
        display_usage
        exit 1
    fi
    DOWNSTREAM_STATS="true"
    CMD+=" --tstat \$SLURM_CPUS_PER_TASK"
elif [[ "$6" == "false" || "$BAM_2" == "-" ]]; then
    CMD+=" --statoff"
    DOWNSTREAM_STATS="false"
else
    printf "Please set a valid value for 6th arg.\n"
    display_usage
    exit 1
fi 

#Paired samples
if [[ -n "$7" && "$7" == "true" ]]; then
    if [[ $DOWNSTREAM_STATS == "false" ]]; then
        printf "Statistical analysis is set to false, therefore the 7th arg (paired analysis) is useless\n"
    fi
    CMD+=" --paired-stats"
fi

#From previous run
if [[ -z "$8" || "$8" == "false" || "$8" == "-" ]]; then
    previous_run="false"
    CMD+=" --task both --tmp \$PWD"
elif [[ "$8" == "post" ]]; then
    #rmats file should be in the output dir already
    CMD+=" --task post --tmp $OUTDIR"
    previous_run="false"
elif [[ "$8" == "maser" ]]; then
    previous_run="true"
else
    printf "Please set a valid value for the 8th argument.\n"
    display_usage
    exit 1
fi

#Statistical thresholds for maser
if [[ -z "$9" || "$9" == "-" ]]; then
    min_avg_reads=5
    fdr_threshold=0.05
    deltaPSI_threshold=0.2
else
    IFS=","
    read -r -a array <<< "${9}"
    min_avg_reads=${array[0]}
    fdr_threshold=${array[1]}
    deltaPSI_threshold=${array[2]}
fi

#Read type
if [[ -z "${10}" || "${10}" == "paired" || "${10}" == "-" ]]; then
    READTYPE="paired"
elif [[ "${10}" == "single" ]]; then
    READTYPE="single"
else
    printf "Please set a valid value for 10th arg.\n"
    display_usage
    exit 1
fi

#Strandness
if [[ -z "${11}" || "${11}" == "fr-firststrand" || "${11}" == "-" ]]; then
    LIBTYPE="fr-firststrand"
elif [[ "${11}" == "fr-secondstrand" ]]; then
    LIBTYPE="$7"
elif [[ "${11}" == "fr-unstranded" ]]; then
    LIBTYPE="${11}"
else
    printf "Please set a valid value for 11th arg.\n"
    display_usage
    exit 1
fi

#Read length
if [[ -z "${11}" || "${11}" == "-" ]]; then
    CMD+=" --readLength 150 --variable-read-length"
else
    CMD+=" --readLength ${11} --variable-read-length"
fi

#--enable-unicode=ucs4
if [[ $READTYPE == "paired" ]];then
    CMD+=" --libType $LIBTYPE"
fi


cat > rmats.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=rmats
#SBATCH --time=72:00:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --image=mcfonsecalab/rmats:latest
#SBATCH --output=%j_rmats.log

echo "------------"
echo "--- CMD ----"
echo "------------"
echo "$CMD"

if [[ $previous_run == "false" ]]; then
    workdir="/home/pedro.barbosa/scratch/rna_seq/splicing/rMATS/\$SLURM_JOB_ID"
    mkdir \$workdir && cd \$workdir
    srun shifter $CMD
fi

if [[ $DOWNSTREAM_STATS == "true" ]]; then
    if [[ $previous_run == "true" ]]; then
        cd $OUTDIR
    fi
    printf "Producing tables of significant events\n"
    printf "Tresholds used for the filtering (fdr, dPSI, min_avg_reads): $fdr_threshold, $deltaPSI_threshold, $min_avg_reads\n"
    echo "srun shifter Rscript /python_env/run_maser.R $label1 $label2 $min_avg_reads $fdr_threshold $deltaPSI_threshold "JC" $GTF"
    #srun shifter Rscript /python_env/run_maser.R $label1 $label2 $min_avg_reads $fdr_threshold $deltaPSI_threshold "JC" $GTF
    #srun cut -f3 coverage_filt/* | sort | uniq | grep -v "geneSymbol" | sed 's/\\(ENSG[0-9]*\\)\\.[0-9]*/\\1/g' > ${OUT_BASENAME}_rmats_negative_list_to_GO.txt
    #srun cut -f2 sign_events_* | sed 's/\\(ENSG[0-9]*\\)\\.[0-9]*/\\1/g' | sort | uniq > ${OUT_BASENAME}_rmats_positive_list_to_GO.txt
    #srun echo "gene_id\tFDR\tdPSI\n" > ${OUT_BASENAME}_rmats_gene_ranks_to_GSEA.txt
    #srun awk '{print \$3, \$(NF-3), \$NF}' OFS="\t" coverage_filt/* |  sed 's/\\(ENSG[0-9]*\\)\\.[0-9]*/\\1/g' | grep -v geneSymbol >> ${OUT_BASENAME}_rmats_gene_ranks_to_GSEA.txt
    #rename 's/sign/${OUT_BASENAME}_sign/g' *
    
    printf "Done\nGenerating sashimi plots from all significant events..\n"
    events=("SE" "A5SS" "A3SS" "RI" "MXE")
    for e in "\${events[@]}"; do
        printf "Looking at significant \${e} events..\n"
        awk -F'\t' 'NR==FNR{c[\$1]++;next};c[\$1]' *sign_events_\${e}.tsv \${e}.MATS.JC.txt > \${e}_toSashimi.txt
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


mv * $OUTDIR
cd ../ && rm -rf \$SLURM_JOB_ID
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
EOL

sbatch rmats.sbatch
