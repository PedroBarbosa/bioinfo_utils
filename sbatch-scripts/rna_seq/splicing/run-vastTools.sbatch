#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=75G
#SBATCH --job-name="vast" 
#SBATCH --output=%j_vasttools.out
#SBATCH --image=biocorecrg/vast-tools:2.1.3
#SBATCH --array=0-0%1


fastq=$(readlink -f "$1")
DEL=$2
p2=${DEL/1/2}
out_dir=$(readlink -f "$3")
if [[ ! -d "${out_dir}/expr_out" &&  ! -d "${out_dir}/to_combine" ]]; then
    mkdir "${out_dir}/expr_out"
    mkdir "${out_dir}/to_combine"
fi
readarray -t pair1 < <(cat $fastq | grep -e "R1.fastq" -e "R1.fq" -e "_1.fastq" -e "_1.fq")

scratch_out="/home/pedro.barbosa/scratch/rna_seq/splicing/$SLURM_JOB_ID"
mkdir $scratch_out
cd $scratch_out

#OUT_BASENAME=$(basename ${pair1[$SLURM_ARRAY_TASK_ID]} | cut -f1,2,3 -d "_")_out.txt
db="/home/rluis/Rui-testing/Genome/hg19_hg38_vast-tools2"
#db="/home/pedro.barbosa/mcfonseca/shared/genomes/mouse/mm9/"
dir=$(dirname ${pair1[$SLURM_ARRAY_TASK_ID]})
pair2=$(basename ${pair1[$SLURM_ARRAY_TASK_ID]/$DEL/$p2})
fullpathpair2="$dir/$pair2"

CMD="vast-tools align ${pair1[$SLURM_ARRAY_TASK_ID]} $fullpathpair2 --dbDir $db --expr -c $CPUS_PER_TASK --sp Hsa" # --output $OUT_BASENAME"

echo $CMD
srun shifter $CMD

cd vast_out
mv expr_out/* ${out_dir}/expr_out
mv to_combine/* ${out_dir}/to_combine
#cd ../../ && rm -rf $SLURM_JOB_ID

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID
