#!/bin/bash
display_usage(){
 printf "Script to automatically generate a sbatch script with a specific header to run simple commands within lobo.\n
 Usage:
    -1st argument must be sbatch script name.
    -2nd argument must be the sbatch job name.
    -3rd argument must be the memory available for the job (in Gb)
    -4th argument is optional. Refers to the number of nodes,tasks and cpus per task, respectively, to employ on this slurm job in lobo (tasks will be set in parallel,not in the srun command). Default:1,1,10. '-' skips this argument.\n"
}

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
    printf "Error, please set the required arguments for the script.\n"
    display_usage
    exit
fi

re='^[0-9]+$'
##JOB SETTINGS""
if [[ -z "$4" || "$9" = "-" ]]; then
    NODES=1
    NTASKS=1
    CPUS_PER_TASK=10
else
    IFS=','
    read -r -a array <<< "$4"
    if [ ${#array[@]} = 3 ]; then
        for elem in "${array[@]}"
        do
            if ! [[ "$elem" =~ $re ]]; then
                printf "Error. Please set INT numbers for the number of nodes, tasks and cpus per task.\n"
                display_usage
                exit 1
            fi
        done
        NODES=${array[0]}
        NTASKS=${array[1]}
        CPUS_PER_TASK=${array[2]}
    else
        printf "ERROR. 3 fields are required for the 4th argument (nodes,tasks,cpus per task). You set a different number.\n"
        display_usage
        exit 1
    fi
fi


cat > /home/pedro.barbosa/$1 <<EOL
#!/bin/bash
#SBATCH --job-name=$2
#SBATCH --time=24:00:00
#SBATCH --mem=${3}G
#SBATCH --nodes=$NODES
#SBATCH --ntasks=$NTASKS
#SBATCH --cpus-per-task=$CPUS_PER_TASK
#SBATCH --workdir=/home/pedro.barbosa/scratch
#SBATCH --output=/home/pedro.barbosa/scratch/%j_${2}.log
EOL
