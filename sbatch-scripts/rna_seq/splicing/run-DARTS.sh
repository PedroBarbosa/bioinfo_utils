#!/usr/bin/env bash
display_usage() {
echo 'Script to run DARTS based on the output from rMATS.

-1st argument must be the directory where rmats output results are stored.
-2nd argument must be the output directory.
-3rd argument is optional. Whether informative prior file should be generated for Darts BHT procedure. If set to true, Darts_DNN model will be employed to predict the probability of differential splicing based on exon-specific cis-features and sample-specific trans features. Default: false. Values: [true|false|-]'
}


if [ -z "$1" ] || [ -z "$2" ]; then
    printf "ERROR:Please provide at least the 2 first arguments required for the script.\n\n"
    display_usage 
    exit 1
fi

DIR=$(readlink -f "$1")
OUT_DIR=$(readlink -f "$2")
if [[ ! -d "$OUT_DIR" ]]; then
    mkdir $(readlink -f "$OUT_DIR")
fi

if [[ -z "$3" || "$3" == "false" ]]; then
    CMD="Darts_BHT bayes_infer -c 0.2 --nthread 10 --od $OUT_DIR --verbose 1"
fi


cat > darts.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=darts
#SBATCH --time=72:00:00
#SBATCH --mem=240G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --image=mcfonsecalab/rmats:latest
#SBATCH --output=%j_darts.log

#scratch_dir="/home/pedro.barbosa/scratch/rna_seq/splicing/darts/\${SLURM_JOB_ID}"
#mkdir \$scratch_dir && cd \$scratch_dir
cd $DIR
events=("SE" "A5SS" "A3SS" "RI" "MXE")
for e in "\${events[@]}"; do
    printf "Looking at significant \${e} events..\n"
    FINAL_CMD="$CMD --rmats-count JC.raw.input.\${e}.txt --annot fromGTF.\${e}.txt -t \${e}"
    printf "\$FINAL_CMD"
    srun shifter \$FINAL_CMD
done

mv * $OUT_DIR
cd ../ && rm -rf \$SLURM_JOB_ID
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
EOL

sbatch darts.sbatch
