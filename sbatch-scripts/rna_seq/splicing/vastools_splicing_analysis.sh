#!/bin/bash
display_usage(){
 printf "Script to automatically run different vastools utilities after aligning the reads.
Merged and non-merged analysis will be performed. In the first, we are more stringent in the compare filtering, while in the non-merged we are more flexible (--noVLOW is removed and --min_range = 5, instead of 10).\n\n
 Usage:
    -1st argument must the directory where the parent directory of the 'to_combine' is located. (e.g. output dir of the align script)
    -2nd argument must be the groups file to perform comparison and/or merging. 1st columns must be the sample name, 2nd column the group it belongs.
    -3rd argument must be the name of the comparison to use in the output files.
    -4th argument is optional. Referes whether the comparison is paired or not. Default:false. Values:[true|false|-]
    -5th argument is optional. Refers to the vastDB to be used. Default: Human hg38. Set '-' to ignore this argument. 
\n"
    
}

if [ -z "$1" ] | [ -z "$2" ] | [ -z "$3" ]; then
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
NAME="$3"

if [[ -z "$4" || "$4" == "-" || "$4" == "false" ]]; then
    PAIRED="false"
elif [[ "$4" == "true" ]]; then
    PAIRED="true"
else
    printf "Please set a valid value for the 4th argument.\n"
    display_usage
    exit 1
fi

if [[ -z "$5" || "$5" == "-" ]]; then
    #vastDB="/home/rluis/Rui-testing/Genome/hg19_hg38_vast-tools2"
    #species="Hsa"
    #assembly="hg38"
    vastDB="/home/mcfonseca/shared/genomes/mouse/mm9/"
    species="Mmu"
    assembly="mm9"
else
    vastDB=$(readlink -f "$5")
    species=""
    assembly""
    printf "You set a specific database. Which species does that refer ?\n"
    display_usage
    exit 1
fi


cat > vast_splicing_analysis.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=vast_splicing
#SBATCH --time=72:00:00
#SBATCH --mem=75G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=%j_vast_splicing.log
#SBATCH --image=biocorecrg/vast-tools:2.1.3

cd $DIR

if [[ ! -d "expr_out" &&  ! -d "to_combine" ]]; then
    printf "Required directory with output of vast align does not exist in the directory ($DIR) you provided.\n"
    display_usage
    exit 1
fi

############
###MERGE####
############
mkdir merged_${NAME}
mkdir merged_${NAME}/to_combine
CMD_merge="vast-tools merge -g $GROUP --sp $species --dbDir $vastDB -o $DIR"
printf "\n###CMD merge##: \${CMD_merge}.\n"
srun shifter \$CMD_merge

declare -a GROUPS_ARRAY
prev_group=""
samples_first_group=""
colors_to_plot=("blue" "red")
i=1
printf "Order\tSampleName\tGroupName\tRColorCode\n" > groups_to_plot_${NAME}.txt
printf "Order\tSampleName\tGroupName\tRColorCode\n" > merged_${NAME}/groups_to_plot_${NAME}.txt
while read line;do
    sample=\$(echo \${line} | awk '{print \$1}')
    group=\$(echo \$line | awk '{print \$2}')
    if [[ \${samples_first_group} == "" ]]; then 
        samples_first_group="\${sample}"
        GROUP_ARRAY+=( "\$group" )
        prev_group="\$group"
        second_group="false"
        printf "\$i\t\$sample\t\$group\t\${colors_to_plot[0]}\n" >> groups_to_plot_${NAME}.txt
        printf "1\t\$group\t\$group\t\${colors_to_plot[0]}\n" >> merged_${NAME}/groups_to_plot_${NAME}.txt
    elif [[ "\$group" != "\$prev_group" ]]; then
        GROUP_ARRAY+=( "\$group" )
        samples_second_group="\${sample}"
        prev_group="\$group"
        second_group="true"
        printf "\$i\t\$sample\t\$group\t\${colors_to_plot[1]}\n" >> groups_to_plot_${NAME}.txt
        printf "2\t\$group\t\$group\t\${colors_to_plot[1]}\n" >> merged_${NAME}/groups_to_plot_${NAME}.txt
    elif [[ \$second_group == "false" ]];then 
        samples_first_group="\${samples_first_group},\${sample}"
        printf "\$i\t\$sample\t\$group\t\${colors_to_plot[0]}\n" >> groups_to_plot_${NAME}.txt
    else
        samples_second_group="\${samples_second_group},\${sample}"  
        printf "\$i\t\$sample\t\$group\t\${colors_to_plot[1]}\n" >> groups_to_plot_${NAME}.txt
    fi
    i=\$((\$i+1))
done < $GROUP

mv to_combine/\${GROUP_ARRAY[0]}\.* merged_${NAME}/to_combine
mv to_combine/\${GROUP_ARRAY[1]}\.* merged_${NAME}/to_combine


###########
##COMBINE##
###########
#No merged
CMD_combine="vast-tools combine -sp $species Hsa -a $assembly --dbDir $vastDB -o $DIR"
printf "\n##CMD combine:## \${CMD_combine}\n"
srun shifter \$CMD_combine
INCLUSION_FILE="0_${NAME}_INCLUSION_LEVELS.tab"
mv INCLUSION_LEVELS* \$INCLUSION_FILE

#Merged
printf "\nCombining now the merged samples..\n"
cd merged_${NAME}
CMD_combine="vast-tools combine -sp $species Hsa -a $assembly --dbDir $vastDB -o \$PWD"
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
CMD_tidy="vast-tools tidy \$INCLUSION_FILE_MERGED -min_Fr 0.5 -outFile 1_${NAME}_merged_TIDY.tsv --p_IR --add_names"
printf "\n##CMDcTIDY MERGED:## \${CMD_tidy}\n"
srun shifter \$CMD_tidy

CMD_compare="vast-tools compare \$INCLUSION_FILE_MERGED -a \${GROUP_ARRAY[0]} -b \${GROUP_ARRAY[1]} --sp $species --dbDir $vastDB \
--min_dPSI 20 \
--min_range 10 \
--max_dPSI 5 \
--noVLOW \
--GO \
--p_IR --use_int_reads --fr_int_reads 0.4"
   
if [[ $PAIRED == "true" ]]; then
    CMD_compare="\${CMD_compare} --paired"
fi
CMD_COMPARE_TO_PLOT="srun shifter \${CMD_compare} --outFile 2_${NAME}_compare_noDPSI_toPlot.tsv"
CMD_COMPARE_WITH_DPSI="srun shifter \${CMD_compare} --print_dPSI --outFile 2_${NAME}_compare_withDPSI.tsv"
\$CMD_COMPARE_TO_PLOT
\$CMD_COMPARE_WITH_DPSI
mkdir 2_aux_files && mv AltEx*.txt BG*.txt IR*txt 2_aux_files
printf "\n##COMPARE MERGED:## \$CMD_COMPARE_TO_PLOT\n"


###################
#Non merged groups#
###################
cd ../
CMD_tidy="vast-tools tidy \$INCLUSION_FILE -min_Fr 0.5 -outFile 1_${NAME}_TIDY.tsv --p_IR --add_names"
printf "\n##TIDY NOT MERGED:## \${CMD_tidy}\n"
srun shifter \$CMD_tidy

CMD_compare="vast-tools compare \$INCLUSION_FILE -a \${samples_first_group} -b \${samples_second_group} --sp $species --dbDir $vastDB \
--min_dPSI 20 \
--min_range 5 \
--max_dPSI 5 \
--GO \
--p_IR --use_int_reads --fr_int_reads 0.4"
   
if [[ $PAIRED == "true" ]]; then
    CMD_compare="\${CMD_compare} --paired"
fi
CMD_COMPARE_TO_PLOT="srun shifter \$CMD_compare --only_samples --outFile 2_${NAME}_compare_noDPSI_toPlot.tsv"
CMD_COMPARE_WITH_DPSI="srun shifter \$CMD_compare --print_dPSI --outFile 2_${NAME}_compare_withDPSI.tsv"
\$CMD_COMPARE_TO_PLOT
\$CMD_COMPARE_WITH_DPSI
mkdir 2_${NAME}_aux_files && mv AltEx*.txt BG*.txt IR*txt 2_${NAME}_aux_files
printf "\n##COMPARE NON MERGED GROUPS:## \$CMD_COMPARE_TO_PLOT\n"


#########################
####Plot from compare####
#########################
srun shifter vast-tools plot -u TRUE -v TRUE -c groups_to_plot_${NAME}.txt -m 2000 -o \$PWD 2_${NAME}_compare_noDPSI_toPlot.tsv
srun shifter vast-tools plot -u TRUE -v TRUE -c merged_${NAME}/groups_to_plot_${NAME}.txt -m 2000 -o merged_${NAME} merged_${NAME}/2_${NAME}_compare_noDPSI_toPlot.tsv
mv \${2_${NAME}_compare_noDPSI_toPlot.tsv/tsv/PSI_plots.pdf} 3_${NAME}_PSI_plots_from_compare.pdf 
mv merged_${NAME}/\${2_${NAME}_compare_noDPSI_toPlot.tsv/tsv/PSI_plots.pdf} merged_${NAME}/3_${NAME}_merged_PSI_plots_from_compare.pdf 
echo "\${2_${NAME}_compare_noDPSI_toPlot.tsv/tsv/PSI_plots.pdf} 3_${NAME}_PSI_plots_from_compare.pdf"
echo "merged_${NAME}/\${2_${NAME}_compare_noDPSI_toPlot.tsv/tsv/PSI_plots.pdf} merged_${NAME}/3_${NAME}_merged_PSI_plots_from_compare.pdf"

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
#srun shifter vast-tools diff -a \${samples_first_group} -b \${samples_second_group} -i \$INCLUSION_FILE \
#            --minDiff 0.2 \
#            --minSamples 3 \
#            --prob 0.95 \
#            --minReads 20 \
#            -c 10 
#            -o \$PWD

echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
EOL

sbatch vast_splicing_analysis.sbatch
