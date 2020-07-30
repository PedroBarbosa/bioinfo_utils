#!/bin/bash
display_usage(){
 printf "Script to automatically run majiq PSI for a set of files.\n
 Usage:
    -1nd argument must the majiq files to calculate PSI. Names will be extracted from filenames directly.
    -2rd argument must be the output directory.
    -3th is optional. Refers to the minimum number of experiments in a group to consider an LSV as valid. Default: (0.5, half of samples). Can be used as a exact number of samples (int value) Use '-' to skip argument and apply default.
    -4rd is optional. It refers to the group name. If set, all majiq files as will be treated as replicates (or samples) from the same group. Use '-' to skip argument and apply default. Default: each file is a sample.\n"

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

majiq_files=$(readlink -f "$1")
OUT=$(readlink -f "$2")
if [[ ! -d "$OUT" ]];then
    mkdir $OUT
fi

if [[ -z "$3" || "$3" == "-" ]]; then
    n_experiments=0.5
else
    n_experiments="$3"
fi
 
files=""
CMD="majiq psi --mem-profile --min-experiments $n_experiments"


if [[ -z "$4" ]];then
    grouped="false"
    declare -a RUN_ARRAY
    while read line; do
        sample_name="$(echo $(basename $line .majiq) | cut -f1,2 -d "_")"
        cmd="$CMD --logger ${sample_name}_psi.log --name $sample_name -o . $line"
        RUN_ARRAY+=("$cmd")
   done < $majiq_files
   printf '%s\n' "${RUN_ARRAY[@]}" > $PWD/psi_parallel_cmds.txt
  
else
    grouped="true"
    group_name="$4"
    while read line; do        
        files="$files $line"
    done < $majiq_files
    CMD="$CMD --logger ${group_name}_psi.log --name $group_name -o \$PWD $files"
fi

cat > psi_majiq.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=majiq_psi
#SBATCH --time=72:00:00
#SBATCH --mem=120G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=%j_majiq_psi.log
#SBATCH --image=mcfonsecalab/majiq:latest 

scratch_out=/home/pedro.barbosa/scratch/rna_seq/majiq/\$SLURM_JOB_ID
mkdir \$scratch_out
cd \$scratch_out

if  [ $grouped == "false" ];then
    while read line;do
        srun shifter \$line --nproc \$SLURM_CPUS_PER_TASK
    done < $PWD/psi_parallel_cmds.txt
#    cat $PWD/psi_parallel_cmds.txt | parallel -j \$SLURM_CPUS_PER_TASK "srun shifter --image=mcfonsecalab/majiq:latest {}"
#    printf '%s\n' "${RUN_ARRAY[@]}" #| parallel -j 10 echo {}
#    parallel -j \$SLURM_CPUS_PER_TASK "srun shifter {}" ::: "${RUN_ARRAY[@]}"
else
    printf "##PSI CMD##\n$CMD\n"
    srun shifter $CMD
fi
#mv * $OUT
#cd ../ && rm -rf \$SLURM_JOB_ID
#conda deactivate
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
EOL

sbatch psi_majiq.sbatch
