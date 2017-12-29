#!/bin/bash
display_usage(){
 printf "Script to automatically run fastqc in lobo from a set of raw paired-end fastq files.\n
 Usage:
    -1st argument must the fastq files to process. Pairs must come consecutively in the file.
    -2nd argument must be the output directory.\n"
}

if [ -z "$1" ] || [ -z "$2" ]; then
    printf "Error. Please set the required parameters\n"
    display_usage
    exit 1
fi

workdir="/home/pedro.barbosa/scratch/fastqc"
if [ ! -d $workdir ]; then
    mkdir $workdir
fi

if [ ! -d $2 ]; then
    outdir=$(readlink -f "$2")
    mkdir $outdir
else
    outdir=$(readlink -f "$2")
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


cat > $workdir/runFastqc.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=qc
#SBATCH --output="$workdir/%j_fastqc.out"
#SBATCH --mail-user=pedro.barbosa@medicina.ulisboa.pt
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=1
#SBATCH --workdir=$workdir
#SBATCH --image=argrosso/htspreprocessing:0.1.1
#SBATCH --export=ALL

srun="srun -N1 -n1"
parallel="parallel -j \$SLURM_NTASKS_PER_NODE --delay 0.2 --joblog parallel.log"
cat $fastq_data | \$parallel '\$srun shifter fastqc -o $workdir -t \$SLURM_CPUS_PER_TASK {}'
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID 
mv $workdir/* $outdir
EOL

sbatch $workdir/runFastqc.sbatch 
