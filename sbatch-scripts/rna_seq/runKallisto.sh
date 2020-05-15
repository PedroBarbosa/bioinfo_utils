#!/usr/bin/env bash
display_usage() {
echo 'Script to run Kallisto for multiple fastq files.
Read groups are automatically added to each output based on the sample basename.

-1st argument must be the file listing RNA-seq pairs consecutively. One file per line.
-2nd argument must be the file of the reference indexed database. (Default hg38 gencode v33 with spike-ins)
-3rd argument must be the output directory.
-4th argument is optional. It is the identifier to extract the sample pair names from fastq files. Default: "_1.fastq"'
}


if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] ; then
        printf "ERROR:Please provide at least the 3 first arguments required for the script.\n\n"
        display_usage 
        exit 1
fi

readarray FASTQ < $(readlink -f "$1")
JOBS=$(( ${#FASTQ[@]} / 2 ))
if [[ $2 == "-" ]]; then
    INDEX="/home/mcfonseca/shared/genomes/human/hg38/kallisto/gencode_v33_with_spike_ins/gencode.v33.with.spikes.idx"
else
    INDEX=$(readlink -f "$2")
fi


if [[ ! -d $(readlink -f "$3") ]]; then
    mkdir $(readlink -f "$3")
fi
OUT=$(readlink -f "$3")

if [[ -z "$4" ]]; then
    DEL="_1.fastq"
    p2="_2.fastq"
else
    DEL="$4"
    p2=${DEL/1/2}
fi

cat > kallisto.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=kallisto
#SBATCH --array=0-$(( $JOBS - 1 ))%5
#SBATCH --time=72:00:00
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --image=mcfonsecalab/rna_quantification:latest
#SBATCH --output=%j_kallisto.log

readarray -t pair1 < <(cat \$(readlink -f "$1") | grep -e "R1.fastq" -e "R1.fq" -e "_1.fastq" -e "_1.fq")
scratch_out=/home/pedro.barbosa/scratch/rna_seq/\$SLURM_JOB_ID
mkdir \$scratch_out
cd \$scratch_out
srun="srun -N1 -n1"
OUT_BASENAME=\$(basename \${pair1[\$SLURM_ARRAY_TASK_ID]} | cut -f1,2,3 -d "_")

dir=\$(dirname \${pair1[\$SLURM_ARRAY_TASK_ID]})
pair2=\$(basename \${pair1[\$SLURM_ARRAY_TASK_ID]/$DEL/$p2})
fullpathpair2="\$dir/\$pair2"
HDF5_USE_FILE_LOCKING=FALSE
\$srun shifter kallisto quant -i $INDEX -o \$PWD -b100 -t \$SLURM_CPUS_PER_TASK \${pair1[\$SLURM_ARRAY_TASK_ID]} \$fullpathpair2
mv abundance.tsv \${OUT_BASENAME}_abundance.tsv
mv abundance.h5 \${OUT_BASENAME}.h5
mv *json \${OUT_BASENAME}_run_info.json
mv * $OUT
cd ../ && rm -rf \$SLURM_JOB_ID
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
EOL

sbatch kallisto.sbatch
