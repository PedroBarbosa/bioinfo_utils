#!/usr/bin/env bash
display_usage() {
echo 'Script to run rMATs for multiple BAM files.

-1st argument must be the file listing the replicates for condition 1. Replicates must be split by ",".
-2st argument must be the file listing the replicates for condition 2. Replicates must be split by ","
-3rd argument must be the GTF annotation file (must be unzipped).
-4th argument must be the output directory.
-5th argument must be the labels for each group, split by ",".
-6th argument is optional. Refers to the read type. Values [paired|single|-]. Default: paired.
-7th argument is optional. Refers to the library type of the RNA. Values: [fr-firstsrand|fr-secondstrand|fr-unstranded|-]. Default: fr-firststrand.
-8th argument is optional. Refers to whether statistical analysis should be performed. Values: [true|false|-]. Default: true'

}


if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] ; then
        printf "ERROR:Please provide at least the 4 first arguments required for the script.\n\n"
        display_usage 
        exit 1
fi

BAM_1=$(readlink -f "$1")
BAM_2=$(readlink -f "$2")
GTF=$(readlink -f "$3")
if [[ ! -d $(readlink -f "$4") ]]; then
    mkdir $(readlink -f "$4")
fi
OUTDIR=$(readlink -f "$4")

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

if [[ -z "$8" || "$8" == "true" ]]; then
    DO_STAT_ANALYSIS="true"
elif [[ "$8" == "false" ]]; then
    DO_STAT_ANALYSIS="false"
else
    printf "Please set a valid value for 8th arg.\n"
    display_usage
    exit 1
fi 

CMD="rmats4.py --b1 $BAM_1 --b2 $BAM_2 --gtf $GTF --od \$PWD -t paired --nthread \$SLURM_CPUS_PER_TASK --cstat 0.0001"
#CMD="Darts_BHT rmats_count --b1 $BAM_1 --b2 $BAM_2 --gtf $GTF --od \$PWD -t paired --nthread \$SLURM_CPUS_PER_TASK"

#--enable-unicode=ucs4
if [[ $LIBTYPE == "paired" ]];then
    CMD="$CMD -libType $LIBTYPE"
fi

if [[ $DO_STAT_ANALYSIS == "false" ]]; then
    CMD="$CMD --statoff"
fi

cat > rmats.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=rmats
#SBATCH --time=72:00:00
#SBATCH --mem=240G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --image=mcfonsecalab/rmats:latest
#SBATCH --output=%j_rmats.log

workdir="/home/pedro.barbosa/scratch/rna_seq/splicing/\$SLURM_JOB_ID"
mkdir \$workdir && cd \$workdir
srun shifter $CMD

printf "Producing tables of significant events\n"
srun shifter Rscript /python_env/run_maser.R $label1 $label2 

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

mv * $OUTDIR
cd ../ && rm -rf \$SLURM_JOB_ID
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
EOL

sbatch rmats.sbatch
