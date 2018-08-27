#!/bin/bash
display_usage(){
 printf "Script to automatically convert a list of bam files into fastq. It will use samtools for the job.\n
 Usage:
    -1st argument must be the list of BAM files to process.
    -2nd argument must refer to the alignment type: PE, SE. 
    -3nd argument is optional. Referes to the output directory to store all fastq files. [default: output is written on the bam directory.]\n"
}

if [ -z "$1" ] || [ -z "$2" ]; then
    printf "Error. Please set the required parameters\n"
    display_usage
    exit 1
fi
while read line
do
    if [ ! -e "$line" ]; then
        printf "Error. $line doesn't exist. Please check more carefully files passed in the 1st argument.\n"
        display_usage
        exit 1
    fi
done < "$1"
bamfiles=$(readlink -f "$1")
base_workdir="/home/pedro.barbosa/scratch/bamConvert"
if [ ! -d $base_workdir ]; then
    mkdir $base_workdir
fi

if [ "$2" == "PE" ]; then
    isPE=true
elif [ "$2" == "SE" ]; then
    isPE=false
else
    printf "2nd argument must either be SE or PE.\n" 
    exit 1
fi

out_dir=""
if [ -z "$3" ]; then
    out_dir="{//}"
elif [ -d "$3" ]; then
    echo "$3 directory will be used to write the final files"
    out_dir=$(readlink -f "$3")
elif [ ! -d "$3" ]; then
    printf "$3 directory doesn't exist. Please set a valid one!\n"
    display_usage
    exit 1
fi


cat > $base_workdir/bam2fastq.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=bam_convert
#SBATCH --output=/home/pedro.barbosa/scratch/bamConvert/%j_convert.log
#SBATCH --mail-user=pedro.barbosa@medicina.ulisboa.pt
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --mem=200G
#SBATCH --nodes=2
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=5
#SBATCH --workdir=$base_workdir
#SBATCH --image=mcfonsecalab/htstools_plus:latest

workdir="$base_workdir/\$SLURM_JOB_ID"
mkdir \$workdir
cd \$workdir

srun="srun --exclusive -N1 -n1 -c5"
parallel="parallel --delay 0.2 -j \$SLURM_NTASKS --joblog parallel.log"

if [ $isPE == "true" ]; then
    cat $bamfiles | \$parallel 'echo "$(date): analysis started!" > {= s{.*/}{};s/\_[^_]+$//;s/\_[^_]+$//; =}.log && \$srun shifter samtools sort -n -@\$SLURM_CPUS_PER_TASK {} | shifter samtools fastq -N -@\$SLURM_CPUS_PER_TASK -s {=s{.*/}{};s/\_[^_]+$//;s/\_[^_]+$//;=}_singletons.fq -1 {= s{.*/}{};s/\_[^_]+$//;s/\_[^_]+$//; =}_R1.fq -2 {= s{.*/}{};s/\_[^_]+$//;s/\_[^_]+$//; =}_R2.fq - && echo "$(date): fastq files generated." >> {= s{.*/}{};s/\_[^_]+$//;s/\_[^_]+$//; =}.log && gzip {= s{.*/}{};s/\_[^_]+$//;s/\_[^_]+$//; =}_R1.fq {= s{.*/}{};s/\_[^_]+$//;s/\_[^_]+$//; =}_R2.fq {=s{.*/}{};s/\_[^_]+$//;s/\_[^_]+$//;=}_singletons.fq && echo "$(date): fastq files gzipped. All done!!" >> {= s{.*/}{};s/\_[^_]+$//;s/\_[^_]+$//; =}.log && mv {= s{.*/}{};s/\_[^_]+$//;s/\_[^_]+$//; =}*gz $out_dir'
else
    cat $bamfiles | \$parallel 'echo "$(date): analysis started!" > {= s{.*/}{};s/\_[^_]+$//;s/\_[^_]+$//; =}.log && \$srun shifter samtools sort -n -@\$SLURM_CPUS_PER_TASK {} | shifter samtools fastq -@\$SLURM_CPUS_PER_TASK -0 {=s{.*/}{};s/\_[^_]+$//;s/\_[^_]+$//;=}.fq  - && echo "$(date): fastq file generated." >> {= s{.*/}{};s/\_[^_]+$//;s/\_[^_]+$//; =}.log && gzip  {=s{.*/}{};s/\_[^_]+$//;s/\_[^_]+$//;=}.fq  && echo "$(date): fastq file gzipped. All done!!" >> {= s{.*/}{};s/\_[^_]+$//;s/\_[^_]+$//; =}.log && mv {=s{.*/}{};s/\_[^_]+$//;s/\_[^_]+$//;=}.fq.gz $out_dir'

fi
    
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID 

EOL

sbatch $base_workdir/bam2fastq.sbatch
