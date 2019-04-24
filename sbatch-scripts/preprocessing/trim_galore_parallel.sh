#!/bin/bash
display_usage(){
 printf "Script to automatically run trim galore in lobo for a set of raw paired-end fastq files.\n
 Usage:
    -1st argument must the fastq files to process. Pairs must come consecutively in the file.
    -2nd argument must be the minimum Phred quality score to trim low quality ends from reads (in addition to adpaters).
    -3rd argument must be the minimum length allowed in a read after its processing.
    -4th argument must be set to turn on/off FastQC on the fastq file once trimming is complete.Values [true|false]. 
    -5th argument must be the output directory. If '-' is set, output directory will be the same as the parent of the input fastq
    -6th argument is optional. Refers to the stringency value. Default: 5.
    -7th argument is optional. If set, turns off/on GNU parallel. Values [true|false]. Default: true.
    -8th argument is optional. Refers to the adapter sequence to trim. By default, trim_galore auto detects it. Values: [DNA sequence| -]
    -9th argument is optional. Refers to the number of bp to remove from the 5' end of the reads. Values: [int| -]. '-' skips the argument.
    -10th argument is optional. Refers to the number of bp to remove from the 3' end of the reads after the adaptor removal [unwanted bias at 3' not
related with adaptors]\n"
}

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ]; then
    printf "Error. Please set the required parameters\n"
    display_usage
    exit 1
fi

workdir="/home/pedro.barbosa/scratch/trim_galore"
cmd="trim_galore"
re="^[0-9]+$"
if [ ! -z "$9" -a "$9" != "-" ]; then
    if ! [[ $9 =~ $re ]] ; then
        echo "error: $9 not a number" 
        display_usage
        exit 1
    else
        cmd="$cmd --clip_R1 "$9" --clip_R2 "$9""
    fi
fi

if [ ! -z "${10}" -a "${10}" != "-" ]; then
    if ! [[ ${10} =~ $re ]] ; then
        echo "error: ${10} not a number" 
        display_usage
        exit 1
    else
        cmd="$cmd --three_prime_clip_R1 "${10}" --three_prime_clip_R2 "${10}""
    fi
fi

if [ ! -z "$8" -a "$8" != "-" ]; then
    adaptor="-a $8" 
    cmd="$cmd $adaptor" 
fi

if [ -z "$7" ] || [ "$7" = "true" ]; then
    runParallel=true
elif [ "$7" = "false" ];then
    runParallel=false
else
    printf "Please set a valid value for the 6th argument."
    display_usage
    exit 1
fi
if [ -z "$6" ]; then
    cmd="$cmd --stringency 5"
else
    cmd="$cmd --stringency $6"
fi

if [ "$4" = "true" ];then
    cmd="$cmd --fastqc"
elif [ "$4" != "false" ]; then
    echo "Please set a valid value for the 3rd argument"
    exit 1
fi

if [ ! -d $workdir ]; then
    mkdir $workdir
fi

if [ "$5" == "-" ]; then
    outdir="{//}"
elif [ ! -d $5 ]; then
    outdir=$(readlink -f "$5")
    mkdir $outdir
else
    outdir=$(readlink -f "$5")
fi

while read -r file; do
   if [[ -e "$file" ]]; then
      continue 
   else
      echo "$file does not exist."
      exit 1
   fi
done < "$1"
fastq_data=$(readlink -f "$1")
cmd="$cmd --gzip --length $3 --trim-n --max_n 2 --paired -q $2"


cat > $workdir/runTrimGalore.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=trimGalore
#SBATCH --mail-user=pedro.barbosa@medicina.ulisboa.pt
#SBATCH --time=72:00:00
#SBATCH --mem=100G
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=6
##SBATCH --workdir=$workdir
#SBATCH --image=docker:argrosso/htspreprocessing:0.1.1
#SBATCH --export=ALL
#SBATCH --output=%j_trimGalore.log


scratch_workdir="$workdir/\$SLURM_JOB_ID"
mkdir \$scratch_workdir
cd \$scratch_workdir

echo \$scratch_workdir
if [ "$runParallel" = "true" ]; then
	srun="srun --exclusive -N1 -n1"
	parallel="parallel --delay 0.2 -j \$SLURM_NTASKS --joblog parallel.log --resume-failed"
        time cat "$fastq_data" | grep -v "R2" | \$parallel 'mkdir {/.}; cd {/.}; \$srun shifter $cmd {} {=s/R1/R2/=}; mv * $outdir; cd ../; rm -rf {/.}'
else
	for file1 in \$(find \$(dirname $fastq_data) -type f -name \$(basename $fastq_data) -exec cat {} \; | grep -v "R2"); 
		do f2=\$(echo \$file1 | sed 's/R1/R2/');
		time srun -N1 -n1 --exclusive shifter $cmd \$file1 \$f2; done
                if [ $outdir = {//} ]; then
                    mv \${file1}* \${f2}* \$(dirname \$file1)
                fi
fi
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID 
if [ $outdir != {//} ]; then
    mv * $outdir
fi
EOL
sbatch $workdir/runTrimGalore.sbatch 
