#!/bin/bash
display_usage(){
 printf "Script to run strelka2 pipeline in lobo for a bunch of bam files.\n
 Usage:.
    -1st argument must be the bam the tumor alignment files. A '.bai' index for each bam file must exist in each bam directory.
    -2nd argument must be the bam normal alignment files. If "-" is set, normal alignment file will be skipped.
    -3nd argument must be the name of the final ouptut directory.
    -4rd argument must refer to the type of analysis in hand. Values:[WGS|exome|targeted].
    -5th argument must be the fasta reference for which the reads were aligned. If '-' is set, the existing hg38 version in lobo will be used. Be aware that a fasta index file (.fai) must also be present in the fasta directory.
    -6th argument is optional. It refers to a bgzip-compressed and tabix-indexed bed file that contains the regions to restrict variant calling. Required when 'exome' or 'targeted' are passed in the 4th argument. If you don't have it yeat, use bgzip and tabix software to create such files.\n"
}

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] ||  [ -z "$4" ]; then
    printf "Error. Please set the required parameters of the script\n"
    display_usage
    exit 1
fi

####CHECK BAM INPUT####
if [ ! -f "$1" ]; then
    echo "File "$1" not found!"
    exit 1
fi
TUMOR_DATA=$(readlink -f "$1")

normalSkip="False"
if [ -f "$2" ]; then
    NORMAL=$(readlink -f "$2")
elif [ "$2" != "-" ]; then
    echo "Please set a valid value for the 2nd argument"
    exit 1
else
    normalSkip="True"
fi
####SCRATCH WORKDIR####
OUTDIR=$(readlink -f "$3")
WORKDIR="/home/pedro.barbosa/scratch/strelka2"
if [ ! -d $WORKDIR ];then
    mkdir $WORKDIR
fi
if [ ! -d $OUTDIR ];then
    mkdir $OUTDIR
fi

###MODE####
CMD="configureStrelkaSomaticWorkflow.py"
analysis=("WGS" "exome" "targeted")
if [[ !  " ${analysis[@]} " =~ " ${4} " ]]; then
    printf "Please set a valid value for the 4th argument.\n"
    display_usage
    exit 1
elif [ "$4" != "WGS" ]; then
    CMD="$CMD --${4}"
    if [ -z "$6" ]; then
        printf "When exome or targeted sequencing, you should provide the target regions in the 5th argument as a single bed bgzip-compressed and tabix-indexed file.\n"
        display_usage
        exit 1
    elif [[ ${6: -4} == ".bgz" && -f "${6}.tbi" ]]; then
        CMD="$CMD --callRegions="$(readlink -f $6)""
    else
        printf "Invalid bed input for target regions. Please make sure you provide a bgzipped bed file (extension .bgz) in the 5th argument, and that in the same directory exists a tabix index of the bgz file.\n"
        display_usage
        exit 1
    fi  
fi

##REFERENCE##
if [ "$5" = "-" ]; then
    CMD="$CMD --reference=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/WholeGenomeFasta/genome.fa"
elif [ ! -f "$5" ]; then
    printf "Please provide a valid fasta file in th 5th argument.\n"
    display_usage
    exit 1
elif [ ! -f "${5}.fai" ]; then
    printf "Fasta index ${5}.fai not found in the reference directory. Please create one with samtools faidx.\n"
    display_usage
    exit 
else
    CMD="$CMD --reference=$5 --outputCallableRegions"
fi

CMD_ROOT=$CMD
##CMD EXTENSION##
if [ "$normalSkip" == "True" ]; then
    CMD="$CMD --tumorBam $TUMOR_DATA"
else
    CMD="$CMD --tumorBam $TUMOR_DATA --normalBam $NORMAL"
fi

echo "$CMD"
cat > $WORKDIR/configureStrelka.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=strelka2
#SBATCH --time=72:00:00
#SBATCH --mem=240G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --image=mcfonsecalab/strelka:2.8.4
#SBATCH --workdir=$WORKDIR
#SBATCH --output=$WORKDIR/%j_strelka2.log

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
\$srun shifter $CMD
echo "\$(timestamp) -> Configuration completed. Will run now workflow."
runWorkflow="\$scratch_out/StrelkaSomaticWorkflow/runWorkflow.py -m local -j \$SLURM_CPUS_ON_NODE"
\$srun shifter \$runWorkflow
echo -e "\$(timestamp) -> Finished job."
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
echo -e "\$(timestamp) -> All done!"
cd ../ && mv \$SLURM_JOB_ID \$SLURM_JOB_ID_strelka2.log $OUTDIR
EOL

sbatch $WORKDIR/configureStrelka.sbatch
sleep 1 
cd $WORKDIR
mv configureStrelka.sbatch $(ls -td -- */ | head -n 1) 

