#!/usr/bin/env bash
display_usage() {
echo 'Script to run rMATs for multiple BAM files.

-1st argument must be the file listing the replicates for condition 1. Replicates must be split by ",".
-2st argument must be the file listing the replicates for condition 2. Replicates must be split by ","
-3rd argument must be the GTF annotation file.
-4th argument must be the output directory.
-5th argument is optional. Refers to the read type. Values [paired|single|-]. Default: paired.
-6th argument is optional. Refers to the library type of the RNA. Values: [fr-firstsrand|fr-secondstrand|fr-unstranded|-]. Default: fr-firststrand.'
}


if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] ; then
        printf "ERROR:Please provide at least the 4 first arguments required for the script.\n\n"
        display_usage 
        exit 1
fi

BAM_1=$(readlink -f "$1")
BAM_2=$(readlink -f "$2")
GTF=$(readlink -f "$3")
OUTDIR=$(readlink -f "$4")

if [[ -z "$5" || "$5" == "paired" || "$5" == "-" ]]; then
    READTYPE="paired"
elif [[ "$5" == "single" ]]; then
    READTYPE="single"
else
    printf "Please set a valid value for 5th arg.\n"
    display_usage
    exit 1
fi

if [[ -z "$6" || "$6" == "fr-firststrand" || "$6" == "-" ]]; then
    LIBTYPE="fr-firststrand"
elif [[ "$6" == "fr-secondstrand" ]]; then
    LIBTYPE="$6"
elif [[ "$6" == "fr-unstranded" ]]; then
    LIBTYPE="$6"
else
    printf "Please set a valid value for 6th arg.\n"
    display_usage
    exit 1
fi

CMD="rmats.py --b1 $BAM_1 --b2 $BAM_2 --gtf $GTF --od \$PWD -t paired --nthread \$SLURM_CPUS_PER_TASK --cstat 0.0001 --nthread 10"

#--enable-unicode=ucs4
if [[ $LIBTYPE == "paired" ]];then
    CMD="$CMD -libType $LIBTYPE"
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
mv * $OUTDIR
cd ../ && rm -rf \$SLURM_JOB_ID
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
EOL

sbatch rmats.sbatch
