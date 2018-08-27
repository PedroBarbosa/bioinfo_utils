#!/bin/bash
display_usage(){
 printf "Script to automatically calculate genome/exome/targeted data coverage.\n
 IMPORTANT: Sample names will be extracted up to the first '_' found on each bam name, make sure that substring is unique.\n
 Usage:.
    -1st argument must be the file containing the bam sequences to analyze. All samples must have been aligned against the same reference.[one per line].
    -2nd argument must be the name of the ouptut directory.
    -3rd argument must be the reference fasta sequence that was used in the alignments. If '-' is set, the hg38 reference present in lobo will be used.
    -4th argument must refer to the type of analysis in hand. Values:[WGS|targeted].
    -5th argument must be provided when the data comes from a targeted experiment. Refers to the target regions in bed format. Use '-' to skip this argument.
    -6th argument is optional. Refers to the bait regions used for the hybridization capture method. If you don't have such file, the targets file will be used as the baits file as well. Use '-' to skip this argument.
    -7th argument is optional. Refers to the maximum distance between a read and the nearest probe/bait/amplicon for the read to be considered 'near probe' and included in percent selected. Default: 250.
    -8th argument is optional. Refers to the max coverage limit for theoritical sensitivity calculations. Default: 600.
    -9th argument is optional. Refers to the number of nodes,tasks and cpus per task, respectively, to employ on this slurm job in lobo (tasks will be set in parallel,not in the srun command). Default:1,5,8. '-' skips this argument.\n"
}

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] ; then
    printf "Error. Please set the required parameters of the script\n"
    display_usage
    exit 1
fi
###PLOTING EXEC
plottingExec="/home/pedro.barbosa/bioinfo_utils/python-scripts/lobo/coverageAnalysis/plotCoverageFromGATK4output.py"

####CHECK BAM INPUT####
if [ ! -f "$1" ]; then
    printf "Error. '$1' file doesn't exist\n"
    exit 1 
    display_usage
fi
while read line
do
    if [ ! -e "$line" ]; then
         printf "Error. $line doesn't exist. Please check more carefully files passed in the 1st argument.\n"
         display_usage
         exit 1
    fi
done < "$1"
BAM_DATA=$(readlink -f "$1")

####SCRATCH WORKDIR####
OUTDIR=$(readlink -f "$2")
WORKDIR="/home/pedro.barbosa/scratch/covAnalysis"
if [ ! -d $WORKDIR ];then
    mkdir $WORKDIR
fi
if [ ! -d $OUTDIR ];then
    mkdir $OUTDIR
fi

##REFERENCE##
if [ "$3" = "-" ]; then
    reference="/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/WholeGenomeFasta/genome.fa"
elif [ ! -f "$3" ]; then
    printf "Please provide a valid fasta file in th 4th argument.\n"
    display_usage
    exit 1
elif [ ! -f "${3}.fai" ]; then
    printf "Fasta index ${3}.fai not found in the reference directory. Please create one with samtools faidx.\n"
    display_usage
    exit 
else
    reference=$(readlink -f "$3")
fi

##NEAR DISTANCE##
re='^[0-9]+$'
if [[ -z "$7" || "$7" = "-" ]]; then
    near_dist="250"
elif ! [[ "$7" =~ $re ]]; then
    printf "Near distance is not a INT number. Please set one valid one in the 7th argument.\n"
    exit 1
    display_usage
else
    near_dist="$7"
fi

##COV CAP##
if [[ -z "$8" || "$8" = "-" ]]; then
    coverage_cap="600"
elif ! [[ "$8" =~ $re ]]; then
    printf "Coverage cap is not a INT number. Please set one valid one in the 8th argument.\n"
    exit 1
    display_usage
else
    coverage_cap="$8"
fi

##JOB SETTINGS""
if [[ -z "$9" || "$9" = "-" ]]; then
    NODES=1
    NTASKS=8
    CPUS_PER_TASK=5
else
    IFS=','
    read -r -a array <<< "$9"
    if [ ${#array[@]} = 3 ]; then
        for elem in "${array[@]}"
        do
            if ! [[ "$elem" =~ $re ]]; then 
                printf "Error. Please set INT numbers for the number of nodes, tasks and cpus per task.\n"
                display_usage
                exit 1
            fi
        done
        NODES=${array[0]}
        NTASKS=${array[1]}
        CPUS_PER_TASK=${array[2]}
    else 
        printf "ERROR. 3 fields are required for the 9th argument (nodes,tasks,cpus per task). You set a different number.\n"
        display_usage
        exit 1
    fi
fi

###MODE####
analysis=("WGS" "targeted")
if [[ !  " ${analysis[@]} " =~ " ${4} " ]]; then
    printf "Please set a valid value for the 4th argument.\n"
    display_usage
    exit 1
elif [ "$4" = "targeted" ]; then
    if [ -z "$5" ]; then
        printf "When targeted sequencing, you should provide the target regions in the 5th argument, and if you have a bed file of the capture baits, you may also provide it in the 6th argument. Please set at least the 5th argument with the target regions.)\n"
        display_usage
        exit 1
    elif [[ ${5: -4} == ".bed" ]]; then
        targets=$(readlink -f "$5")
        if [[ -z "$6" || "$6" = "-" ]] ; then
            baits=$targets
        elif [[ ! -z "$6" && ${6: -4} == ".bed" ]]; then
            baits=$(readlink -f "$6")
        else
            printf "If provided, 6th argument must be in bed format. Please correct this issue.\n"
            exit 1
            display_usage
        fi
        cat > $WORKDIR/targetCoverage.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=target_coverage
#SBATCH --time=72:00:00
#SBATCH --mem=245G
#SBATCH --nodes=$NODES
#SBATCH --ntasks=$NTASKS
#SBATCH --cpus-per-task=$CPUS_PER_TASK
#SBATCH --image=broadinstitute/gatk:latest
#SBATCH --workdir=$WORKDIR
#SBATCH --output=$WORKDIR/%j_getTargetCoverage.log

ulimit -c ulimited
timestamp() {
    date +"%Y-%m-%d  %T"
}
export -f timestamp
scratch_out=$WORKDIR/\$SLURM_JOB_ID
mkdir \$scratch_out
cd \$scratch_out
srun="srun -N1 -n1 --slurmd-debug 2"
parallel="parallel --delay 0.2 -j \$SLURM_NTASKS --env timestamp --joblog parallel.log --resume-failed"
echo "\$(timestamp) -> Analysis started! Converting bed files to Picard IntervalList format."
seq_dict=\$(head -n 1 $BAM_DATA)
echo "\$(timestamp) -> \$seq_dict file used to extract sequence dictionary required to GATK BedToIntervalList utility"
tg=\$(basename ${targets})
bt=\$(basename ${baits})
\$srun shifter gatk BedToIntervalList -I=$targets -O=\${tg/.bed/_t.picard} -SD=\$seq_dict
\$srun shifter gatk BedToIntervalList -I=$baits -O=\${bt/.bed/_b.picard} -SD=\$seq_dict

echo "\$(timestamp) -> Calculating individual coverages and experiment metrics!"
header_hsmetrics="sample\tbait_region\ttarget_region\t#reads_in_bam\tfraction_locatedOnOrNearBaits\tfraction_away_baits\tfraction_alignedBP_offTarget\tmean_cov_targets\tmedian_cov_targets\tfraction_targetsNoCov\tfraction_targetBP_1X_cov\tfraction_2X_cov\tfraction_10X_cov\tfraction_20X_cov\tfraction_30X_cov\tfraction_40X_cov\tfraction_50X_cov\tfraction_100X_cov"
echo -e "\$header_hsmetrics" > final_collectHSmetrics_all.txt

BAITS=\${bt/.bed/_b.picard}
TARGETS=\${tg/.bed/_t.picard}
cmd="-BI=\$BAITS -TI=\$TARGETS --MINIMUM_BASE_QUALITY=15 --MINIMUM_MAPPING_QUALITY=10 --METRIC_ACCUMULATION_LEVEL=ALL_READS --COVERAGE_CAP=$coverage_cap --NEAR_DISTANCE=$near_dist -R=$reference"
##CMD="gatk CollectHsMetrics -BI=\$BAITS -TI=\$TARGETS --interval_merging=OVERLAPPING_ONLY --MINIMUM_BASE_QUALITY=15 --MINIMUM_MAPPING_QUALITY=10 --METRIC_ACCUMULATION_LEVEL=ALL_READS --COVERAGE_CAP=$coverage_cap --NEAR_DISTANCE=$near_dist -R=$reference"

##PARALLEL
echo \$srun
echo \$BAITS
cat $BAM_DATA | \$parallel 'i=\$(basename {}); out=\$(echo \$i | cut -f1 -d "_"); \$srun shifter gatk CollectHsMetrics -BI final_gencode.v19_HCM_concatenated_noRepeatsOnIntrons_b.picard -TI final_gencode.v19_HCM_concatenated_noRepeatsOnIntrons_t.picard --MINIMUM_BASE_QUALITY=15 --MINIMUM_MAPPING_QUALITY=10 --METRIC_ACCUMULATION_LEVEL=ALL_READS --COVERAGE_CAP=$coverage_cap --NEAR_DISTANCE=$near_dist -R=$reference --INPUT={} --OUTPUT=\${out}_HS_metrics.txt --PER_TARGET_COVERAGE=\${out}_perTargetCov.txt' 


#awk '/BAIT_SET/{getline; print}' \${out}_HS_metrics.txt | awk 'BEGIN{OFS="\t";} {print \$3,\$4,\$6,\$19,\$20,\$34,\$23,\$24,\$29,\$36,\$37,\$38,\$39,\$40,\$41,\$42,\$43}' >> final_collectHSmetrics_all.txt; sed -i '\$s/^/'"\${out}\t"'/' final_collectHSmetrics_all.txt' 

#unique_samples=()
#for j in \$(find $BAM_DATA -exec cat {} \; );do
#    printf "\$(timestamp): Processing \$(basename \$j) file.\n"
#    i=\$(basename \$j)
#    out=\$(echo \$i | cut -f1 -d "_")
#    if [[ " \${unique_samples[@]} " =~ " \${out} " ]]; then
#        printf "\$(timestamp): Warning, duplicate sample ID (\${out}) found after splitting by first '_'. Will use the filename instead without the extension."
#        out="\${i%.*}"
#    fi
#    unique_samples+=(\${out})
#    \$srun shifter \$CMD -I=\$j -O=\${out}_HS_metrics.txt --PER_TARGET_COVERAGE=\${out}_perTargetCov.txt
#    awk '/BAIT_SET/{getline; print}' \${out}_HS_metrics.txt | awk 'BEGIN{OFS="\t";} {print \$3,\$4,\$6,\$19,\$20,\$34,\$23,\$24,\$29,\$36,\$37,\$38,\$39,\$40,\$41,\$42,\$43}' >> final_collectHSmetrics_all.txt
#    sed -i '\$s/^/'"\${out}\t"'/' final_collectHSmetrics_all.txt

#    echo -e "##\${out}" >> final_perTargetCoverage.txt
#    cat \${out}_perTargetCov.txt >> final_perTargetCoverage.txt
#    printf "\$(timestamp): Done!\n"
#done

if [ -f "$plottingExec" ];then
    printf "\$(timestamp): Plotting some of the results!"
    \$srun shifter --image=mcfonsecalab/python36_bio:latest python $plottingExec \$scratch_out targeted
else
    printf "\$(timestamp): No plotting script found."
fi
##CalculateTargetCoverage was removed from last GATK4 release. Piece of code removed.
#echo -e "\\n\$(timestamp) -> Calculating proportional read counts per target and sample!"
#CMD="gatk CalculateTargetCoverage -L \$TARGETS -O final_pCovCounts.txt --transform PCOV -targetInfo FULL -groupBy SAMPLE --rowSummaryOutput final_perFeature_totalCounts.txt --columnSummaryOutput final_perSample_totalTargetsCoverage.txt"
#while read line; do
#    CMD="\$CMD -I \$line"
#done < $BAM_DATA
#\$srun shifter \$CMD
#mv * ../\$SLURM_JOB_ID*log $OUTDIR
EOL
        sbatch $WORKDIR/targetCoverage.sbatch
        sleep 1
        cd $WORKDIR
        mv targetCoverage.sbatch $(ls -td -- */ | head -n 1)
    else
        printf "Invalid input for target regions. Please provide a bed file.\n"
        display_usage
        exit 1
    fi  
elif [ "$4" = "WGS" ]; then
    if [[ -z "$5" || "$5" = "-" ]] && [[ -z "$6" || "$6" = "-" ]];then
        printf "Whole genome coverage analysis will be employed.\n"
    elif [[ ${5: -4} == ".bed" && -f $5 ]] && [[ -z "$6" || "$6" = "-" ]]; then
        printf "Coverage analysis will be restricted to the regions provided in the "$5" file.\n"
        targets=$(readlink -f "$5")
    else
        printf "Error. Either you didn't provide a correct bed file in the 5th argument or you wrongly set the 6th argument, which shouldn't be passed on WGS experiments (don't set it or skip with '-'.\n"
        display_usage
        exit 1
    fi
    cat > $WORKDIR/wgsCoverage.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=wgs_coverage
#SBATCH --time=72:00:00
#SBATCH --mem=245G
#SBATCH --nodes=$NODES
#SBATCH --ntasks=$NTASKS
#SBATCH --cpus-per-task=$CPUS_PER_TASK
#SBATCH --image=broadinstitute/gatk:latest
#SBATCH --workdir=$WORKDIR
#SBATCH --output=$WORKDIR/%j_getWGScoverage.log

ulimit -c ulimited
timestamp() {
    date +"%Y-%m-%d  %T"
}
export -f timestamp
scratch_out=$WORKDIR/\$SLURM_JOB_ID
mkdir \$scratch_out
cd \$scratch_out
srun="srun -N1 -n1 --slurmd-debug 2"
parallel="parallel -k --delay 0.2 -j \$SLURM_NTASKS --env timestamp --joblog parallel.log --resume-failed"

if [ -n "$targets" ]; then   
    echo "\$(timestamp) -> Analysis started! Since a bed file was provided to restrict the coverage assessment, a conversion of such file into a Picard IntervalList format needs to be performed."
    seq_dict=\$(head -n 1 $BAM_DATA)
    echo "\$(timestamp) -> \$seq_dict file used to extract sequence dictionary required to GATK BedToIntervalList utility"
    tg=\$(basename ${targets})
    \$srun shifter gatk BedToIntervalList -I=$targets -O=\${tg/.bed/.picard} -SD=\$seq_dict
    INTERVALS="--INTERVALS=\${tg/.bed/.picard}"
fi
echo "\$(timestamp) -> Calculating individual coverages and experiment metrics!"
header_wgsmetrics="sample\tgenome_territory\tmean_coverage\tstd_coverage\tmedian_coverage\tfraction_alignedBp_excluded\tfraction_refBp_1X_cov\tfraction_refBp_5X_cov\tfraction_refBp_10X_cov\tfraction_refBp_20X_cov\tfraction_refBp_30X_cov\tfraction_refBp_40X_cov\tfraction_refBp_50X_cov\tfraction_refBp_70X_cov\tfraction_refBp_90X_cov\tfraction_refBp_100_cov"
echo -e "\$header_wgsmetrics" > final_colllectWgsMetrics_all.txt

CMD="gatk --java-options '-Xmx245G' CollectWgsMetrics --MINIMUM_BASE_QUALITY=15 --MINIMUM_MAPPING_QUALITY=10 --COVERAGE_CAP=$coverage_cap --INCLUDE_BQ_HISTOGRAM=true -R=$reference"
#Again, parallel with issues
#cat "$BAM_DATA" | \$parallel '\$srun shifter -I={} -O={/.}_WGS_metrics.txt \$CMD \$INTERVALS'
unique_samples=()
for j in \$(find $BAM_DATA -exec cat {} \; );do
    printf "\$(timestamp): Processing \$(basename \$j) file.\n"
    i=\$(basename \$j)
    out=\$(echo \$i | cut -f1 -d "_")
    if [[ " \${unique_samples[@]} " =~ " \${out} " ]]; then
        printf "\$(timestamp): Warning, duplicate sample ID (\${out}) found after splitting by first '_'. Will use the filename instead without the extension."
        out="\${i%.*}"
    fi
    unique_samples+=(\${out})
    \$srun shifter \$CMD -I=\$j -O=\${out}_WGS_metrics.txt \$INTERVALS
    awk '/GENOME_TERRITORY/{getline; print}' \${out}_WGS_metrics.txt | cut -f1,2,3,4,12,13,14,15,17,19,20,21,23,25,26 >> final_colllectWgsMetrics_all.txt
    sed -i '\$s/^/'"\${out}\t"'/' final_colllectWgsMetrics_all.txt
 
    sed -n -e '/^coverage/,\$p' \${out}_WGS_metrics.txt > \${out}_wgscov.histo
    printf "\$(timestamp): Sample \$i processed!\n"
done
echo "\$(timestamp) -> Done!!!"

if [ -f "$plottingExec" ];then
    printf "\$(timestamp): Plotting some of the results!"
    \$srun shifter --image=mcfonsecalab/python36_bio:latest python $plottingExec \$scratch_out WGS
else
    printf "\$(timestamp): No plotting script found."
fi
mv * ../\$SLURM_JOB_ID*log $OUTDIR
EOL
    sbatch $WORKDIR/wgsCoverage.sbatch
    sleep 1
    cd $WORKDIR
    mv wgsCoverage.sbatch $(ls -td -- */ | head -n 1)

fi
