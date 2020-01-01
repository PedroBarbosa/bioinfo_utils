#!/bin/sh
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=120G
#SBATCH --job-name="vasttools" 
#SBATCH --output=slurm_%j_vasttools.out
#SBATCH --mail-user=rluis@medicina.ulisboa.pt
#SBATCH --mail-type=ALL
#SBATCH --image=biocorecrg/vast-tools:2.1.3


CommandToRUN="vast-tools combine -sp Hsa -a hg38 --dbDir /home/rluis/Rui-testing/Genome/hg19_hg38_vast-tools2 -o $1"
echo COMMAND USED - \" $CommandToRUN  \"
srun shifter $CommandToRUN

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID
