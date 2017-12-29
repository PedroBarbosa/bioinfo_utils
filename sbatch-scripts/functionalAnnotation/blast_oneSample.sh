#!/bin/bash
display_usage(){
 printf "Script to automatically run blast in lobo for one single fasta.\n
 Usage:
    -1st argument must be BLAST exec to run.
    -2st argument must be the file containing the sequences to blast.
    -3nd argument must be the path for the database files (except the suffix).
    -4th argument must be the name of the ouptut file.
    -5th argument is optional. If set, refer to a set of parameters for this specific blast run [theads, outfmt and maxTargetSeqs are set by default]\n"
}

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ]; then
    printf "Error. Please set the required parameters\n"
    display_usage
    exit 1
fi

workdir="/home/pedro.barbosa/scratch/blast"

if [ ! -d $workdir ]; then
    mkdir $workdir
fi

cat > $workdir/runBlast.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=blastAlu
#SBATCH --output=/home/pedro.barbosa/scratch/blast/blast-%j.out
#SBATCH --mail-user=pedro.barbosa@medicina.ulisboa.pt
#SBATCH --mail-type=ALL
#SBATCH --time=40:00:00
#SBATCH --mem=70G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --workdir=/home/pedro.barbosa/scratch/blast/
#SBATCH --image=miguelpmachado/blast:2.6.0-broadwell
#SBATCH --export=ALL

time srun shifter $1 -query $2 -db $3 -out /home/pedro.barbosa/scratch/blast/$4 -evalue 1e-5 -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' -num_threads 10 -max_target_seqs 5

echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID 
EOL

sbatch $workdir/runBlast.sbatch
