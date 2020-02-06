#!/bin/bash
display_usage(){
 printf "Script to automatically run bbtools to perform quality control analysis for a set of fastq files.\n
 Usage:
    -1nd argument must the file describing fastq data.
    -2th argument must the file the output directory.
    -3rd argument must be the ribosomal sequences to perform the k-mer base removal. If not given (-), hg38 ribosomal sequences will be used.
    -4th argument is optional. It is the identifier to extract the sample pair names from fastq files. Default: '_R1.fq'.\n"	
}

if [ -z "$1" ] || [ -z "$2" ]; then
    printf "Error. Please set the required parameters\n"
    display_usage
    exit 1
fi

readarray FASTQ < $(readlink -f "$1")
JOBS=$(( ${#FASTQ[@]} / 2 ))

if [[ ! -d $(readlink -f "$2") ]]; then
    mkdir $(readlink -f "$2")
fi
OUT=$(readlink -f "$2")


if [[ -z "$3" || "$3" == "-" ]]; then
    rRNA="/home/pedro.barbosa/mcfonseca/shared/genomes/human/hg38/rRNA_hg38_rseqc.fasta"
else
    rRNA=$(readlink -f "$3")
fi


if [[ -z "$4" ]]; then
    DEL="_R1.fq"
    p2="_R2.fq"
else
    DEL="$4"
    p2=${DEL/1/2}
fi

cat > bbduk.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=bbduk
#SBATCH --time=72:00:00
#SBATCH --mem=120G
#SBATCH --array=0-$JOBS%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=bryce911/bbtools:latest
#SBATCH --output=%j_bbduk.log

readarray -t pair1 < <(cat \$(readlink -f "$1") | grep -e "R1.fastq" -e "R1.fq" -e "_1.fastq" -e "_1.fq")
scratch_out=/home/pedro.barbosa/scratch/rna_seq/\$SLURM_JOB_ID
mkdir \$scratch_out
cd \$scratch_out
srun="srun -N1 -n1"
OUT_BASENAME=\$(basename \${pair1[\$SLURM_ARRAY_TASK_ID]} | cut -f1,2,3 -d "_")

dir=\$(dirname \${pair1[\$SLURM_ARRAY_TASK_ID]})
pair2=\$(basename \${pair1[\$SLURM_ARRAY_TASK_ID]/$DEL/$p2})
fullpathpair2="\$dir/\$pair2"

\$srun shifter bbduk.sh in1=\${pair1[\$SLURM_ARRAY_TASK_ID]} in2=\$fullpathpair2 out1=\${OUT_BASENAME}_bbduk_R1.fq.gz out2=\${OUT_BASENAME}_bbduk_R2.fq.gz ref=$rRNA k=31 
mv * $OUT
cd ../ && rm -rf \$SLURM_JOB_ID
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
EOL

sbatch bbduk.sbatch
