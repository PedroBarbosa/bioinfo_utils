#!/usr/bin/env bash
display_usage() {
echo 'Script to produce sashimi plots from a rMATs run.

-1st argument must be the file listing the replicates for condition 1. Replicates must be split by ",".
-2st argument must be the file listing the replicates for condition 2. Replicates must be split by ","
-3rd argument must be the events file generated from rmats.
-4th argument must be the event type: (SE, RI, MXE, A5SS, A3SS)
-5th argument must be the output directory.
-6th argument must be the labels for each group, split by ",".
'
}


if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ]; then
        printf "ERROR:Please provide at least the 4 first arguments required for the script.\n\n"
        display_usage 
        exit 1
fi

BAM_1=$(cat $(readlink -f "$1"))
BAM_2=$(cat $(readlink -f "$2"))
events=$(readlink -f "$3")
event_type="$4"

if [[ ! -d $(readlink -f "$5") ]]; then
    mkdir $(readlink -f "$5")
fi
OUTDIR=$(readlink -f "$5")

IFS=','
read -r -a array <<< "$6"
label1=${array[0]} 
label2=${array[1]}


CMD="rmats2sashimiplot --b1 $BAM_1 --b2 $BAM_2 -o $OUTDIR --l1 $label1 --l2 $label2 --exon_s 1 --intron_s 5 -t $event_type -e $events"

cat > sashimi_rmats.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=sashimi
#SBATCH --time=72:00:00
#SBATCH --mem=240G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --image=mcfonsecalab/rmats:latest
#SBATCH --output=%j_sashimi.log

workdir="/home/pedro.barbosa/scratch/rna_seq/splicing/\$SLURM_JOB_ID"
mkdir \$workdir && cd \$workdir
echo $CMD
srun shifter $CMD

#while read line; do
#    torun="$CMD -c $line:$GTF"    
#    echo \$torun
#    \$torun
#done < $coordinates_file

#torun="$CMD -c chr10:+:76324400:7655730:$GTF"
#echo \$torun
#srun \$torun


mv * $OUTDIR
cd ../ && rm -rf \$SLURM_JOB_ID
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
EOL

sbatch sashimi_rmats.sbatch
