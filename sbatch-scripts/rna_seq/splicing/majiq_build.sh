#!/bin/bash
display_usage(){
 printf "Script to automatically set up a majiq confi file and run build procedure.\n
 Usage:
    -1nd argument must the file that describes the input data (2cols, 1st filepath, 2nd sample group). Bam location is assumed to be the same for all files. Dirname will be extracted to create the config file.
    -2th argument must the file the output directory.
    -3rd argument is optional. Refers to the gff file. If not given (-), gencode hg38 v34 primary assembly will be used.
    -4th argument is optional. Refers whether to run simplifier to remove non-relevant splicing variations. Default: true. Values:[true|false|-].
    -5th argument is optional. Refers to the genome version to use in the config. Default: 'hg38'. Use '-' to skip argument.
    -6th argument is optional. Refers to the strandness of the RNAseq protocol: Default: 'reverse'. Values:[reverse|forward|None]. Use '-' to skip the argument.
    -7th argument is optional. Refers to the read length of the experiments. Default: 150. Use '-' to skip the argument.
    -8th argument is optional. Refers to the min number of experiments from one group to consider a LSV as valid. Default: 2. Use '-' to skip the argument.
    -9th argument is optional. Refers to the min number of reads combining all positions in a LSV to consider the LSV as existing in data. Default: 5. Use '-' to skip the argument.\n"    
}

if [ -z "$1" ] || [ -z "$2" ]; then
    printf "Error. Please set the required parameters\n"
    display_usage
    exit 1
fi

bams=$(readlink -f "$1")
bam_dirname=$(dirname $(head -1 $bams))
bam_dirname=$(echo $bam_dirname | awk '{print $1}')

#Process samples and groups
list_groups=($(cut -f2 $bams | sort | uniq))
array_string=""
for group in "${list_groups[@]}"; do
    array_string="${array_string}[$group]='' "
done    
declare -A array
arrat=( ${array_string})
while read line; do
    bam_file=$(echo $line | awk '{print $1}')
    group=$(echo $line | awk '{print $2'})
    array[$group]+="$(basename $bam_file .bam),"
done < $bams

#OUTdir
if [[ ! -d $(readlink -f "$2") ]]; then
    mkdir $(readlink -f "$2")
fi
OUT=$(readlink -f "$2")

if [[ -z "$3" || "$3" == "-" ]]; then
    GFF="/home/pedro.barbosa/mcfonseca/shared/genomes/human/hg38/gencode_v34_primary_assembly.gff3"
else
    GFF=$(readlink -f "$3")
fi

CMD="majiq build --mem-profile --logger majiq_build.log"
if [[ -z "$4" || "$4" == "true" || "$4" == "-" ]]; then
    CMD="$CMD --simplify --simplify-annotated 5 --simplify-denovo 3"
fi

if [[ -z "$5" || "$5" == "hg38" || "$5" == "-" ]]; then
    genome="hg38"
else
    genome="$5"
fi

if [[ -z "$6" || "$6" == "reverse" || "$6" == "-" ]]; then
    strandness="reverse"
elif [[ "$6" == "forward" ]]; then
    strandness="forward"
elif [[ "$6" == "None" ]]; then 
    strandness="None"
else
    printf "Please set a valid value for the 6th arg.\n"
    display_usage
    exit 1
fi

if [[ -z "$7" || "$7" == "150" || "$7" == "-" ]]; then
    readlength="150"
else
    readlength="$7"
fi


if [[ -z "$8" || "$8" == "-" || "$8" == 2 ]];then
    minExperiments=2
else
    minExperiments="$8"
fi


if [[ -z "$9" || "$9" == "-" || "$9" == 5 ]];then
    minReads=5
else
    minReads="$9"
fi


cat > majiq_config.txt <<EOL
[info]
readlen=$readlength
bamdirs=$bam_dirname
genome=$genome
strandness=$strandness
[experiments]
EOL
for i in "${!array[@]}"
do
  echo "${i}=${array[$i]::-1}" >> majiq_config.txt
done
config=$(echo $PWD/majiq_config.txt)

cat > build_majiq.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=majiq_build
#SBATCH --time=72:00:00
#SBATCH --mem=120G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --output=%j_majiq_build.log
#SBATCH --image=mcfonsecalab/majiq:latest

scratch_out=/home/pedro.barbosa/scratch/rna_seq/majiq/\$SLURM_JOB_ID
mkdir \$scratch_out
cd \$scratch_out
#source activate majiq
CMD="$CMD -j \$SLURM_CPUS_PER_TASK  --minreads $minReads --min-experiments $minExperiments --min-denovo 7 --min-intronic-cov 0.05 --dump-constitutive --dump-coverage --output \$PWD --conf $config $GFF"
echo \$CMD
srun shifter \$CMD
mv * $OUT
#conda deactivate
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
EOL

sbatch build_majiq.sbatch

