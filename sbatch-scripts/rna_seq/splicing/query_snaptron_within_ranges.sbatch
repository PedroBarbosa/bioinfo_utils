#!/bin/bash
#SBATCH --job-name=query_snaptron
#SBATCH --time=72:00:00
#SBATCH --mem=30Gb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=%j_query_snaptron.log

SNAPTRON_R_SCRIPT="/home/pedro.barbosa/git_repos/bioinfo_utils/r-scripts/rnaseq/splicing/query_snaptron.R"
intervals=$(readlink -f $1)

echo "Processing $intervals file.."
source activate r_env
Rscript $SNAPTRON_R_SCRIPT $intervals
conda deactivate

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID
