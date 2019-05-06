#!/bin/bash
display_usage(){
 printf "Script to run strelka2 pipeline in lobo for a bunch of bam files.\n
 Usage:.
    -1st argument must be the file containing the alignment files. A '.bai' index for each bam file must exist in each bam directory.
    -2nd argument must be the name of the final ouptut directory.
    -3rd argument must refer to the type of analysis in hand. Values:[WGS|exome|targeted].
    -4th argument must be the fasta reference for which the reads were aligned. If '-' is set, the existing hg38 version in lobo will be used. Be aware that a fasta index file (.fai) must also be present in the fasta directory.
    -5th argument is optional. It refers to a bgzip-compressed and tabix-indexed bed file that contains the regions to restrict variant calling. Required when 'exome' or 'targeted' are passed in the 3rd argument. If you don't have it yeat, use bgzip and tabix software to create such files.\n"
}

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
    printf "Error. Please set the required parameters of the script\n"
    display_usage
    exit 1
fi

####CHECK BAM INPUT####
if [ ! -f "$1" ]; then
    echo "File "$1" not found!"
    exit 1
fi
while read line
do
    if [ ! -f "$line" ]; then
         printf "Error. $line doesn't exist. Please check more carefully files passed in the 1st argument.\n"
         display_usage
         exit 1
    fi
done < "$1"
BAM_DATA=$(readlink -f "$1")

####SCRATCH WORKDIR####
OUTDIR=$(readlink -f "$2")
WORKDIR="/home/pedro.barbosa/scratch/strelka2"
if [ ! -d $WORKDIR ];then
    mkdir $WORKDIR
fi
if [ ! -d $OUTDIR ];then
    mkdir $OUTDIR
fi

###MODE####
CMD="configureStrelkaGermlineWorkflow.py"
analysis=("WGS" "exome" "targeted")
if [[ !  " ${analysis[@]} " =~ " ${3} " ]]; then
    printf "Please set a valid value for the 3th argument.\n"
    display_usage
    exit 1
elif [ "$3" != "WGS" ]; then
    CMD="$CMD --${3}"
    if [ -z "$5" ]; then
        printf "When exome or targeted sequencing, you should provide the target regions in the 5th argument as a single bed bgzip-compressed and tabix-indexed file.\n"
        display_usage
        exit 1
    elif [[ ${5: -4} == ".bgz" && -f "${5}.tbi" ]]; then
        CMD="$CMD --callRegions="$(readlink -f $5)""
    else
        printf "Invalid bed input for target regions. Please make sure you provide a bgzipped bed file (extension .bgz) in the 5th argument, and that in the same directory exists a tabix index of the bgz file.\n"
        display_usage
        exit 1
    fi  
fi

##REFERENCE##
if [ "$4" = "-" ]; then
    CMD="$CMD --reference=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/WholeGenomeFasta/genome.fa"
elif [ ! -f "$4" ]; then
    printf "Please provide a valid fasta file in th 4th argument.\n"
    display_usage
    exit 1
elif [ ! -f "${4}.fai" ]; then
    printf "Fasta index ${4}.fai not found in the reference directory. Please create one with samtools faidx.\n"
    display_usage
    exit 
else
    CMD="$CMD --reference=$4"
fi

CMD_ROOT=$CMD
##CMD EXTENSION##
numb_samples=$(wc -l < $BAM_DATA)
numb_samples_to_split=75
batched=false
if [[ $numb_samples -gt $numb_samples_to_split ]]; then
    printf "Number of samples ($numb_samples) exceeds strelka recommendations for germline joint variant calling ($numb_samples_to_split). Will split analysis in several batches.\n"
    samples_array=()
    batched=true
    while read line;do samples_array+=("--bam=${line}"); done < $BAM_DATA
#    printf "Array elements: ${samples_array[*]}.\n"
    batches=$((${#samples_array[@]} / $numb_samples_to_split + 1))
    printf "Will split analysis in $batches batches.\n"
    END=$batches
    index=0
    declare -A dynamic_CMD=()   
    declare -a orders;
    for ((i=1;i<=END;i++)); do
        if [[ $i -eq $batches ]]; then #slice until the end
            CMD="$CMD_ROOT ${samples_array[@]:$index}"
            dynamic_CMD["job_$i"]=$CMD
            orders+=( "job_$i" )
        else
            CMD="$CMD_ROOT ${samples_array[@]:$index:$numb_samples_to_split}"
            index=$(($index + $numb_samples_to_split))
            dynamic_CMD["job_$i"]=$CMD
            orders+=( "job_$i" )
        fi
    done
    printf "Jobs to run: ${!dynamic_CMD[*]}.\n"
    WORKDIR="$WORKDIR/batched_run"
    if [ ! -d $WORKDIR ];then
        mkdir "$WORKDIR"
    fi

else
    printf "Total number of samples is smaller than 50. Will process them all in one batch.\n"
    while read line;do CMD="$CMD --bam=${line}";done < $BAM_DATA
fi

cat > $WORKDIR/configureStrelka.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=strelka2
#SBATCH --time=72:00:00
#SBATCH --mem=240G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --image=mcfonsecalab/strelka:2.8.4
##SBATCH --workdir=$WORKDIR
#SBATCH --output=%j_strelka2.log

timestamp() {
    date +"%Y-%m-%d  %T"
}
export -f timestamp
scratch_out=$WORKDIR/\$SLURM_JOB_ID
mkdir \$scratch_out
cd \$scratch_out
srun="srun -N1 -n1 --slurmd-debug 3"
parallel="parallel --delay 0.2 -j \$SLURM_NTASKS  --env timestamp --joblog parallel.log --resume-failed"
echo "\$(timestamp) -> Job started! Configuring a new workflow running script.."
if [ "$batched" = "true" ];then
    \$srun shifter \$1
else
    \$srun shifter $CMD
fi
echo "\$(timestamp) -> Configuration completed. Will run now workflow."
runWorkflow="\$scratch_out/StrelkaGermlineWorkflow/runWorkflow.py -m local -j \$SLURM_CPUS_ON_NODE"
\$srun shifter \$runWorkflow
echo -e "\$(timestamp) -> Finished job."
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
echo -e "\$(timestamp) -> All done!"
cd ../ mv \$scratch_out \${SLURM_JOB_ID}_strelka2.log $OUTDIR
EOL

if [ "$batched" = "true" ]; then
    echo "#!/bin/bash" > $WORKDIR/runStrelkaInBatches.sh
    chmod +x $WORKDIR/runStrelkaInBatches.sh
    previous_job=""
    for key in "${orders[@]}"
    do
        if [ $key = "job_1" ];then
            echo -e "$key=\$(sbatch $WORKDIR/configureStrelka.sbatch '${dynamic_CMD[$key]}')\n" >> $WORKDIR/runStrelkaInBatches.sh
        else
            echo -e "$key=\$(sbatch --dependency=afterok:\${${previous_job}##* } --kill-on-invalid-dep=yes $WORKDIR/configureStrelka.sbatch '${dynamic_CMD[$key]}')\n" >> $WORKDIR/runStrelkaInBatches.sh
        fi
        previous_job=$key
    done
    $WORKDIR/runStrelkaInBatches.sh
else
    sbatch $WORKDIR/configureStrelka.sbatch
    sleep 1 
    cd $WORKDIR
    mv configureStrelka.sbatch $(ls -td -- */ | head -n 1) 
fi
