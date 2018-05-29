#!/bin/bash
display_usage(){
 printf "Script to automatically run trimmomatic in lobo for a set of raw paired-end fastq files.\n
 Usage:
    -1st argument must the fastq files to process. Pairs must come consecutively in the file.
    -2nd argument must be the minimum Phred quality score to allow in the read in a sliding window trimming approach.
    -3rd argument must be the minimum length allowed in a read after its processing.
    -4th argument must be the output directory. If '-' is set, output directory will be the same as the parent of the input fastq
    -5th argument is optional. Refers to the number of mismatches allowed in the seed alignment against the adapters sequences. Default: 2. Values: [int|-]. Use '-' to skip argument.
    -6th argument is optional. If set, turns off/on GNU parallel. Values [true|false]. Default: true.
    -7th argument is optional. Refers to the number of bp to remove from the 5' end of the reads. Values: [int| -]. '-' skips the argument.
    -8th argument is optional. Refers to the number of bp to remove from the 3' end of the reads after the adaptor removal [unwanted bias at 3' not related with adaptors]\n"
}

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ]; then
    printf "Error. Please set the required parameters\n"
    display_usage
    exit 1
fi

workdir="/home/pedro.barbosa/scratch/trimmomatic"
adapters="/home/pedro.barbosa/resources/illuminaAdpaters/adapters.fasta"
basecmd="java -jar /NGStools/Trimmomatic-0.36/trimmomatic.jar PE"
re="^[0-9]+$"
if [ ! -z "$7" -a "$7" != "-" ]; then
    if ! [[ $7 =~ $re ]] ; then
        echo "error: $7 not a number" 
        display_usage
        exit 1
    else
        args="$args HEADCROP:$7"
    fi
fi

if [ ! -z "$8" -a "$8" != "-" ]; then
    if ! [[ $8 =~ $re ]] ; then
        echo "error: $8 not a number" 
        display_usage
        exit 1
    else
        args="$args CROP:$8"
    fi
fi

if [ -z "$6" ] || [ "$6" = "true" ]; then
    runParallel=true
elif [ "$6" = "false" ];then
    runParallel=false
else
    printf "Please set a valid value for the 6th argument."
    display_usage
    exit 1
fi
if [ -z "$5" ] || [ "$5" == "-" ]; then
    args="$args ILLUMINACLIP:$adapters:2:30:10:8:true"
else
    args="$args ILLUMINACLIP $adapters:$5:30:10:8:true"
fi

if [ ! -d $workdir ]; then
    mkdir $workdir
fi

if [ "$4" == "-" ]; then
    outdir="{//}"
elif [ ! -d $4 ]; then
    outdir=$(readlink -f "$4")
    mkdir $outdir
else
    outdir=$(readlink -f "$4")
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
basecmd="$basecmd -trimlog trimmomatic.log"
args="$args MINLEN:$3 SLIDINGWINDOW:4:$2 LEADING:15 TRAILING:15"

cat > $workdir/runTrimmomatic.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=parallel_trimmomatic
#SBATCH --mail-user=pedro.barbosa@medicina.ulisboa.pt
#SBATCH --time=72:00:00
#SBATCH --mem=230G
#SBATCH --nodes=2
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=3
#SBATCH --workdir=$workdir
#SBATCH --image=ummidock/trimmomatic:0.36-2
#SBATCH --export=ALL
#SBATCH --output=%j_trimmomatic.log


scratch_workdir="$workdir/\$SLURM_JOB_ID"
mkdir \$scratch_workdir
cd \$scratch_workdir

if [ "$runParallel" = "true" ]; then
	srun="srun --exclusive -N1 -n1"
	parallel="parallel --delay 0.2 -j \$SLURM_NTASKS --joblog parallel.log --resume-failed"

        time cat "$fastq_data" | grep -v "R2" | \$parallel 'mkdir {/.}; cd {/.}; basenamef1={/}; f2={=s/R1/R2/=}; basenamef2=\$(basename \$f2);\$srun shifter $basecmd -threads 3 {} \$f2 \${basenamef1/R1/R1_filt} \${basenamef1/R1/R1_filtUnpaired} \${basenamef2/R2/R2_filt} \${basenamef2/R2/R2_filtUnpaired} $args; mv * $outdir; cd ../; rm -rf {/.}'
else
	for file1 in \$(find \$(dirname $fastq_data) -type f -name \$(basename $fastq_data) -exec cat {} \; | grep -v "R2"); 
		do f2=\$(echo \$file1 | sed 's/R1/R2/');
                fout1=\$(basename \$file1)
                fout2=\$(basename \$f2)
		time srun -N1 -n1 --exclusive shifter $basecmd -threads 40 \$file1 \$f2 \${fout1/R1/R1_filt} \${fout1/R1/R1_filtUnpaired} \${fout2/R2/R2_filt} \${fout2/R2/R2_filtUnpaired} $args; done
                if [ $outdir = {//} ]; then
                    mv \${fout1}* \${fout2}* \$(dirname \$file1)
                fi
fi
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID 
if [ $outdir != {//} ]; then
    mv * $outdir
fi
EOL
sbatch $workdir/runTrimmomatic.sbatch 

