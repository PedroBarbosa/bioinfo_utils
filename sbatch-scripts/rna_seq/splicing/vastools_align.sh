#!/bin/bash
display_usage(){
 printf "Script to automatically run vastools align a set of files.\n
 Usage:
    -1nd argument must the fastq files to align against vastools database.
    -2rd argument must be the output directory.
    -3th argument is optional. Refers to the vastDB to be used. Default: Human hg38. Set '-' to ignore this argument. 
    -4th argument is optional. Refers to the delimiter used to distinguish pairs. Default: _R1.fq. Use '-' to skip the argument.
    -5th argument is optional. Refers to the intron retention pipeline to use. Default:1. Values:[1|2|-]
\n"
    
}

if [ -z "$1" ] || [ -z "$2" ]; then
    printf "Error. Please set the required parameters\n"
    display_usage
    exit 1
fi

if [[ ! -f "$(readlink -f $1)" ]]; then
    printf "Error. Please set a valid file in the 1st arg.\n"
    display_usage
    exit 1
fi

fq_files=$(readlink -f "$1")
readarray FASTQ < $(readlink -f "$1")
JOBS=$(( ${#FASTQ[@]} / 2 ))

OUT=$(readlink -f "$2")
if [[ ! -d "$OUT" ]];then
    mkdir $OUT
fi

if [[ -z "$3" || "$3" == "-" ]]; then
    vastDB="/home/mcfonseca/shared/genomes/human/hg38/vast-tools/"
    species="hg38"
    #vastDB="/home/mcfonseca/shared/genomes/mouse/GRCm38.p6/vast-tools/"    
    #species="mm10"
else
    vastDB=$(readlink -f "$3")
    species=""
    printf "You set a specific database. Which species does that refer ?\n"
    display_usage
    exit 1
fi
 

if [[ -z "$4" || "$4" == "-" ]]; then
    DEL="_R1.fq"
    p2="_R2.fq"
else
    DEL="$4"
    p2=${DEL/1/2}
fi

if [[ -z "$5" || "$5" == "-" || "$5" == "1" ]]; then
    IR_version="1"
elif [[ "$5" == "2" ]]; then
    IR_version="2"
else
    printf "Please set a valid IR version pipeline to use.\n"
    display_usage
fi

cat > vast_align.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=vast_align
#SBATCH --time=72:00:00
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --output=%j_vast_align.log
##SBATCH --image=vastgroup/vast-tools:v2.4.0
#SBATCH --image=vastgroup/vast-tools:latest
#SBATCH --array=0-$(( $JOBS -1 ))%26

scratch_out=/home/pedro.barbosa/scratch/rna_seq/splicing/vasttools/\$SLURM_JOB_ID
mkdir \$scratch_out
cd \$scratch_out

if [[ ! -d "${OUT}/expr_out" &&  ! -d "${OUT}/to_combine" ]]; then
#    mkdir "${OUT}/expr_out"
    mkdir "${OUT}/to_combine"
fi

readarray -t pair1 < <(cat $fq_files | grep -e "R1.fastq" -e "R1.fq" -e "_1.fastq" -e "_1.fq")

dir=\$(dirname \${pair1[\$SLURM_ARRAY_TASK_ID]})
pair2=\$(basename \${pair1[\$SLURM_ARRAY_TASK_ID]/$DEL/$p2})
fullpathpair2="\$dir/\$pair2"

#--expr removed. Error in last image
#CMD="vast-tools align \${pair1[\$SLURM_ARRAY_TASK_ID]} \$fullpathpair2 --dbDir $vastDB --IR_version $IR_version --keep -c \$SLURM_CPUS_PER_TASK --sp $species" 
CMD="/home/pedro.barbosa/git_repos/vast-tools/vast-tools align \${pair1[\$SLURM_ARRAY_TASK_ID]} \$fullpathpair2 --dbDir $vastDB --IR_version $IR_version --keep -c \$SLURM_CPUS_PER_TASK --sp $species"
echo \$CMD
srun \$CMD
#srun shifter \$CMD

cd vast_out
#mv expr_out/* ${OUT}/expr_out
mv to_combine/* ${OUT}/to_combine
cd ../../ && rm -rf \$SLURM_JOB_ID

echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID

EOL

sbatch vast_align.sbatch
