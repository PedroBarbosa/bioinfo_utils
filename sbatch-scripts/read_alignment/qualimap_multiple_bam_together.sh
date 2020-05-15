#!/bin/bash
display_usage(){
 printf "Script to automatically run qualimap quality control analysis for a set of bam files.\n
 Usage:
    -1nd argument must the file that describes the input data.
    -2th argument must the file the output directoryi.
    -3rd argument must refer to the gtf file. If not given (-), gencode hg38 v33 primary annotation will be used.\n"
}

if [ -z "$1" ] || [ -z "$2" ]; then
    printf "Error. Please set the required parameters\n"
    display_usage
    exit 1
fi

bams=$(readlink -f "$1")

if [[ ! -d $(readlink -f "$2") ]]; then
    mkdir $(readlink -f "$2")
fi
OUT=$(readlink -f "$2")

if [[ -z "$3" || "$3" == "-" ]]; then
    GTF="/home/pedro.barbosa/mcfonseca/shared/genomes/human/hg38/gencode.v33.primary_assembly.annotation.gtf"
else
    GTF=$(readling -f "$3")
fi

cat > qualimap.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=qualimap
#SBATCH --time=72:00:00
#SBATCH --mem=120G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --image=pegi3s/qualimap:latest
#SBATCH --output=%j_qualimap.log

scratch_out=/home/pedro.barbosa/scratch/rna_seq/\$SLURM_JOB_ID
mkdir \$scratch_out
cd \$scratch_out

srun shifter qualimap multi-bamqc --java-mem-size=8G -data $bams -gff $GTF --run-bamqc -outdir . 

mv * $OUT
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
EOL

sbatch qualimap.sbatch

