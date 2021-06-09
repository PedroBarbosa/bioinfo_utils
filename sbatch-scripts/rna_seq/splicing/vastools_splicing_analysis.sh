#!/bin/bash
display_usage(){
 printf "Script to automatically run different vastools utilities after aligning the reads.
Merged and non-merged analysis will be performed. In the first, we are more stringent in the compare filtering, while in the non-merged we are more flexible (--noVLOW is removed and --min_range = 5, instead of 10).\n\n
 Usage:
    -1st argument must the directory where the parent directory of the 'to_combine' is located. (e.g. output dir of the align script)
    -2nd argument must be the groups file to perform comparison and/or merging. 1st columns must be the sample name, 2nd column the group it belongs.
    -3rd argument must be the output directory of the analysis.
    -4rd argument must be the name of the comparison to use in the output files.
    -5th argument is optional. Referes whether the comparison is paired or not. Default:false. Values:[true|false|-]
    -6th argument is optional. If set, skips merged analysis. Default:false. Merged analysis is performed by default. Values:[true|false|-]
    -7th argument is optional. If set, skips non-merged analysis. Default: false. Non-merged analysis is performed by default. Values:[true|false|-]
    -8th argument is optional. Refers to the vastDB to be used. Default: Human hg38. Set '-' to ignore this argument. 
    -9th argument is optional. Refers to the intron retention pipeline used in the align stage. Default:1. Values:[1|2|-]
\n"
    
}

if [ -z "$1" ] | [ -z "$2" ] | [ -z "$3" ] | [ -z "$4" ]; then
    printf "Error. Please set the required parameters\n"
    display_usage
    exit 1
fi

if [[ ! -d "$(readlink -f $1)" ]]; then
    printf "Error. Please set a valid directory in the 1st arg.\n"
    display_usage
    exit 1
fi

if [[ ! -f "$(readlink -f "$2")" ]];then
    printf "Error. Please provide a valid file in the 2nd arg.\n"
    display_usage
    exit 1
fi 

DIR=$(readlink -f "$1")
GROUP=$(readlink -f "$2")

if [[ ! -d $(readlink -f "$3") ]]; then
    mkdir $(readlink -f "$3")
fi
OUTDIR=$(readlink -f "$3")
NAME="$4"

if [[ -z "$5" || "$5" == "-" || "$5" == "false" ]]; then
    PAIRED="false"
elif [[ "$5" == "true" ]]; then
    PAIRED="true"
else
    printf "Please set a valid value for the 5th argument.\n"
    display_usage
    exit 1
fi

if [[ -z "$6" || "$6" == "-" || "$6" == "false" ]]; then
    DO_MERGE="true"
elif [[ "$6" == "true" ]]; then
    DO_MERGE="false"
else
    printf "Please set a valid value for the 6th argument.\n"
    display_usage
    exit 1
fi

if [[ -z "$7" || "$7" == "-" || "$7" == "false" ]]; then
    DO_NON_MERGE="true"
elif [[ "$7" == "true" ]]; then
    DO_NON_MERGE="false"
else
    printf "Please set a valid value for the 7th argument.\n"
    display_usage
    exit 1
fi



if [[ -z "$8" || "$8" == "-" ]]; then
    vastDB="/home/mcfonseca/shared/genomes/human/hg38/vast-tools/"
    species="Hs2"
    assembly="hg38"
    #vastDB="/home/mcfonseca/shared/genomes/mouse/GRCm38.p6/vast-tools/"
    #species="mm10"
    #species="Mmu"
    #assembly="mm10"
else
    vastDB=$(readlink -f "$8")
    species=""
    assembly=""
    printf "You set a specific database. Which species does that refer ?\n"
    display_usage
    exit 1
fi

if [[ -z "$9" || "$9" == "-" || "$9" == "1" ]]; then
    IR_version="1"
elif [[ "$9" == "2" ]]; then
    IR_version="2"
else
    printf "Please set a valid IR version pipeline to use (9th arg).\n"
    display_usage
fi

cat > vast_splicing_analysis.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=vast_splicing
#SBATCH --time=72:00:00
#SBATCH --mem=75G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --output=%j_vast_splicing.log
##SBATCH --image=vastgroup/vast-tools:v2.4.0
#SBATCH --image=vastgroup/vast-tools:latest

cd $DIR

if [[ ! -d "expr_out" && ! -d "to_combine" ]]; then
    printf "Required directory with output of vast align does not exist in the directory ($DIR) you provided.\n"
    display_usage
    exit 1
fi

############
###MERGE####
############
if [[ $DO_MERGE == "true" ]]; then
    mkdir ${OUTDIR}/merged_${NAME}
    mkdir ${OUTDIR}/merged_${NAME}/to_combine
    CMD_merge="vast-tools merge -g $GROUP --sp $species --dbDir $vastDB -o $DIR --IR_version $IR_version"
    printf "\n###CMD merge##: \${CMD_merge}.\n"
    srun shifter \$CMD_merge
    printf "Order\tSampleName\tGroupName\tRColorCode\n" > $OUTDIR/merged_${NAME}/groups_to_plot_${NAME}.txt
fi

declare -a GROUPS_ARRAY
prev_group=""
samples_first_group=""
colors_to_plot=("blue" "red")
i=1
if [[ $DO_NON_MERGE == "true" ]]; then
    printf "Order\tSampleName\tGroupName\tRColorCode\n" > $OUTDIR/groups_to_plot_${NAME}.txt
fi

while read line;do
    sample=\$(echo \${line} | awk '{print \$1}')
    group=\$(echo \$line | awk '{print \$2}')
    if [[ \${samples_first_group} == "" ]]; then 
        samples_first_group="\${sample}"
        GROUP_ARRAY+=( "\$group" )
        prev_group="\$group"
        second_group="false"
        if [[ $DO_NON_MERGE == true ]]; then
            printf "\$i\t\$sample\t\$group\t\${colors_to_plot[0]}\n" >> $OUTDIR/groups_to_plot_${NAME}.txt
        fi
        if [[ $DO_MERGE == "true" ]]; then
            printf "1\t\$group\t\$group\t\${colors_to_plot[0]}\n" >> $OUTDIR/merged_${NAME}/groups_to_plot_${NAME}.txt
        fi
    elif [[ "\$group" != "\$prev_group" ]]; then
        GROUP_ARRAY+=( "\$group" )
        samples_second_group="\${sample}"
        prev_group="\$group"
        second_group="true"
        if [[ $DO_NON_MERGE == true ]]; then
            printf "\$i\t\$sample\t\$group\t\${colors_to_plot[1]}\n" >> $OUTDIR/groups_to_plot_${NAME}.txt
        fi
        if [[ $DO_MERGE == "true" ]]; then
            printf "2\t\$group\t\$group\t\${colors_to_plot[1]}\n" >> $OUTDIR/merged_${NAME}/groups_to_plot_${NAME}.txt
        fi
    elif [[ \$second_group == "false" ]];then 
        samples_first_group="\${samples_first_group},\${sample}"
        if [[ $DO_NON_MERGE == true ]]; then
            printf "\$i\t\$sample\t\$group\t\${colors_to_plot[0]}\n" >> $OUTDIR/groups_to_plot_${NAME}.txt
        fi
    else
        samples_second_group="\${samples_second_group},\${sample}"  
        if [[ $DO_NON_MERGE == true ]]; then
             printf "\$i\t\$sample\t\$group\t\${colors_to_plot[1]}\n" >> $OUTDIR/groups_to_plot_${NAME}.txt
        fi
    fi
    i=\$((\$i+1))
done < $GROUP

if [[ $DO_MERGE == "true" ]]; then
    mv to_combine/\${GROUP_ARRAY[0]}\.* $OUTDIR/merged_${NAME}/to_combine
    mv to_combine/\${GROUP_ARRAY[1]}\.* $OUTDIR/merged_${NAME}/to_combine
fi


###########
##COMBINE##
###########
if [[ $DO_NON_MERGE == true ]]; then
     #No merged#
     CMD_combine="vast-tools combine -sp $species --dbDir $vastDB -o $DIR --IR_version $IR_version"
     printf "\n##CMD combine:## \${CMD_combine}\n"
     srun shifter \$CMD_combine
     INCLUSION_FILE="0_${NAME}_INCLUSION_LEVELS.tab"
     mv INCLUSION_LEVELS* ${OUTDIR}/\${INCLUSION_FILE}
fi
if [[ $DO_MERGE == "true" ]]; then
    #Merged
    printf "\nCombining now the merged samples..\n"
    cd $OUTDIR/merged_${NAME}
    CMD_combine="vast-tools combine -sp $species --dbDir $vastDB -o \$PWD"
    printf "\n##CMD combine merged:## \$CMD_combine\n"
    srun shifter \$CMD_combine
    INCLUSION_FILE_MERGED="0_${NAME}_merged_INCLUSION_LEVELS.tab"
    mv INCLUSION_LEVELS* \$INCLUSION_FILE_MERGED


    ##########################
    ####TIDY AND COMPARE #####
    ##########################

    #################
    #####Merged######
    #################
    CMD_tidy="vast-tools tidy \$INCLUSION_FILE_MERGED --min_Fr 0.75 --min_SD 0 -outFile 1_${NAME}_merged_TIDY_to_rank_genes.tsv --add_names --log"
    printf "\n##CMDcTIDY MERGED:## \${CMD_tidy}\n"
    srun shifter \$CMD_tidy
    mv *Tidy.log 1_${NAME}_merged_TIDY.log

    CMD_compare="vast-tools compare \$INCLUSION_FILE_MERGED -a \${GROUP_ARRAY[0]} -b \${GROUP_ARRAY[1]} --sp $species --dbDir $vastDB \
--min_dPSI 20 \
--min_range 10 \
--max_dPSI 5 \
--noVLOW \
--noB3 \
--GO \
--p_IR --use_int_reads --fr_int_reads 0.4"
   
    if [[ $PAIRED == "true" ]]; then
        CMD_compare="\${CMD_compare} --paired"
    fi
    CMD_COMPARE_TO_PLOT="srun shifter \${CMD_compare} --outFile 2_${NAME}_compare_noDPSI_toPlot.tsv"
    CMD_COMPARE_WITH_DPSI="srun shifter \${CMD_compare} --print_dPSI --outFile 2_${NAME}_compare_withDPSI.tsv"
    \$CMD_COMPARE_TO_PLOT
    \$CMD_COMPARE_WITH_DPSI
    mkdir 2_aux_files && mv All* AltEx*.txt BG*.txt IR*txt 2_aux_files
    cd 2_aux_files && rm *toPlot.txt && mv All* 2_positive_geneIDs_to_GO.txt && mv BG* 2_negative_geneIDs_to_GO.txt && mv 2_* ../ && cd ../
    printf "\n##COMPARE MERGED:## \$CMD_COMPARE_TO_PLOT\n"
    cd ../
fi

###################
#Non merged groups#
###################
if [[ $DO_NON_MERGE == true ]]; then
    cd $OUTDIR
    CMD_tidy="vast-tools tidy \$INCLUSION_FILE -min_Fr 0.5 --min_SD 0 -outFile 1_${NAME}_TIDY_to_rank_genes.tsv --add_names --log"
    printf "\n##TIDY NOT MERGED:## \${CMD_tidy}\n"
    srun shifter \$CMD_tidy
    mv *Tidy.log 1_${NAME}_TIDY.log

    CMD_compare="vast-tools compare \$INCLUSION_FILE -a \${samples_first_group} -b \${samples_second_group} --sp $species --dbDir $vastDB \
--min_dPSI 20 \
--min_range 5 \
--max_dPSI 5 \
--noB3 \
--GO \
--p_IR --use_int_reads --fr_int_reads 0.4"

    if [[ $PAIRED == "true" ]]; then
        CMD_compare="\${CMD_compare} --paired"
    fi
    CMD_COMPARE_TO_PLOT="srun shifter \$CMD_compare --only_samples --outFile 2_${NAME}_compare_noDPSI_toPlot.tsv"
    CMD_COMPARE_WITH_DPSI="srun shifter \$CMD_compare --print_dPSI --outFile 2_${NAME}_compare_withDPSI.tsv"
    \$CMD_COMPARE_TO_PLOT
    \$CMD_COMPARE_WITH_DPSI
    mkdir 2_${NAME}_aux_files && mv All* AltEx*.txt BG*.txt IR*txt 2_${NAME}_aux_files
    cd 2_${NAME}_aux_files && rm *toPlot.txt && mv All* 2_positive_geneIDs_to_GO.txt && mv BG* 2_negative_geneIDs_to_GO.txt && mv 2_* ../ && cd ../
    printf "\n##COMPARE NON MERGED GROUPS:## \$CMD_COMPARE_TO_PLOT\n"


    #########################
    ####Plot from compare####
    #########################
    srun shifter vast-tools plot -u TRUE -v TRUE -c groups_to_plot_${NAME}.txt -m 2000 -o \$PWD 2_${NAME}_compare_noDPSI_toPlot.tsv
    mv \${2_${NAME}_compare_noDPSI_toPlot.tsv/tsv/PSI_plots.pdf} 3_${NAME}_PSI_plots_from_compare.pdf 
fi

if [[ $DO_MERGE == "true" ]]; then
    srun shifter vast-tools plot -u TRUE -v TRUE -c merged_${NAME}/groups_to_plot_${NAME}.txt -m 2000 -o merged_${NAME} merged_${NAME}/2_${NAME}_compare_noDPSI_toPlot.tsv
    mv merged_${NAME}/\${2_${NAME}_compare_noDPSI_toPlot.tsv/tsv/PSI_plots.pdf} merged_${NAME}/3_${NAME}_merged_PSI_plots_from_compare.pdf 
fi


########################
####Plot from TIDY #####
########################
#Rscript filter_vastools_tidy.R TIDY_inclusion.tsv
#cut -f2 3_TIDY_INCLUSION_dPSI_0.2.tsv | grep -ve "NA" -ve "event" | sort | uniq > 4_EVENTS_with_dPSI_0.2.txt 

#fgrep -wf 4_EVENTS_with_dPSI_0.2.txt 1_INCLUSION_LEVELS_FULL-Hsa3-hg38.tab > 5_INCLUSION_LEVELS_FULL_dPSI_0.2.tab
#head -1 1_INCLUSION_LEVELS_FULL-Hsa3-hg38.tab | cat - 5_INCLUSION_LEVELS_FULL_dPSI_0.2.tab > 5_INCLUSION_LEVELS_FULL_dPSI_0.2.tab2
#mv 5_INCLUSION_LEVELS_FULL_dPSI_0.2.tab2 5_INCLUSION_LEVELS_FULL_dPSI_0.2.tab
#srun shifter vast-tools plot -u TRUE -v TRUE -c $groups_to_plot -m 2000 -o \$PWD 5_INCLUSION_LEVELS_FULL_dPSI_0.2.tab
#cat DiffAS-Hsa4-hg38-dPSI15-range5-p_IR_t0-vs-t12-with_dPSI.tab | awk -vOFS="\t" 'NF{NF-=1}1;' > DiffAS_to_plot.tab


###################
##### DIFF ########
###################
if [[ $DO_NON_MERGE == true ]]; then
    DIFF_CMD="vast-tools diff -a \${samples_first_group} -b \${samples_second_group} \
--sampleNameA \${GROUP_ARRAY[0]} --sampleNameB \${GROUP_ARRAY[1]} \
--i \${INCLUSION_FILE} \
-o \$PWD \
-d 3_${NAME}_DIFF \
--prob 0.9 \
--minDiff 0.2 \
--minSamples 3 \
-c 10"
# --minReads 20 

    if [[ $PAIRED == "true" ]]; then
        DIFF_CMD="\${DIFF_CMD} --paired=TRUE"
    fi
    printf "\n##DIFF NON MERGED GROUPS:## \$DIFF_CMD\n"
    srun shifter \$DIFF_CMD
    awk '\$6 >= 0.2' 3_${NAME}_DIFF.tab > 3_${NAME}_DIFF_FILT.tab
#    head -1 3_${NAME}_DIFF.tab | cat - tmp > 3_${NAME}_DIFF_FILT.tab && rm tmp
fi

if [[ $DO_MERGE == "true" ]]; then
    cd merged_${NAME}
    DIFF_CMD_MERGED="vast-tools diff -a \${GROUP_ARRAY[0]} -b \${GROUP_ARRAY[1]} \
-i \${INCLUSION_FILE_MERGED} \
-o \$PWD \
-d 3_${NAME}_DIFF_MERGED \
--prob 0.9 \
--minDiff 0.2 
--minReads 20 \
-c 10"

    printf "\n##DIFF MERGED GROUPS:## \$DIFF_CMD_MERGED\n"
    srun shifter \$DIFF_CMD_MERGED
    awk '\$6 >= 0.2' 3_${NAME}_DIFF_MERGED.tab > 3_${NAME}_DIFF_MERGED_FILT.tab
#    head -1 3_${NAME}_DIFF_MERGED.tab | cat - tmp > 3_${NAME}_DIFF_MERGED_FILT.tab && rm tmp
fi

rm -rf raw*
cd ../../ && rm -rf raw*

echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
EOL

sbatch vast_splicing_analysis.sbatch
