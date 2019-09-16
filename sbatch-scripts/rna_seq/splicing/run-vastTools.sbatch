#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=240G
#SBATCH --job-name="vast" 
#SBATCH --output=%j_vasttools.out
#SBATCH --image=nunoagostinho/vast-tools:2.0.2
#SBATCH --array=0-3%4


fastq=$(readlink -f "$1")
DEL=$2
p2=${DEL/1/2}

readarray -t pair1 < <(cat $fastq | grep -e "R1.fastq" -e "R1.fq" -e "_1.fastq" -e "_1.fq")

scratch_out=/home/pedro.barbosa/scratch/rna_seq/splicing
cd $scratch_out

#OUT_BASENAME=$(basename ${pair1[$SLURM_ARRAY_TASK_ID]} | cut -f1,2,3 -d "_")_out.txt

#db="/home/rluis/Rui-testing/Genome/hg19_hg38_vast-tools2"
db="/home/pedro.barbosa/mcfonseca/shared/genomes/mouse/mm9/"
dir=$(dirname ${pair1[$SLURM_ARRAY_TASK_ID]})
pair2=$(basename ${pair1[$SLURM_ARRAY_TASK_ID]/$DEL/$p2})
fullpathpair2="$dir/$pair2"

CMD="vast-tools align ${pair1[$SLURM_ARRAY_TASK_ID]} $fullpathpair2 --dbDir $db --expr -c 20 --sp Mmu" # --output $OUT_BASENAME"

echo $CMD
srun shifter $CMD

#mv $basename $SLURM_SUBMIT_DIR

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID
