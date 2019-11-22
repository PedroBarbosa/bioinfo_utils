#!/bin/bash
display_usage(){
 printf "Script to automatically run RNA seq quality control analysis for a set of bam files.\n
 Usage:
    -1st argument must the bam files to process.
    -2nd argument must the output directory.
    -3rd argument is optional. Refers to the annotation in bed format. Default(hg38 gencode v31 obtained with bedops). Use '-' to skip this argument.\n"
}

if [ -z "$1" ] || [ -z "$2" ]; then
    printf "Error. Please set the required parameters\n"
    display_usage
    exit 1
fi

bams=$(readlink -f "$1")
JOBS=$(cat $bams | wc -l)
if [[ -z "$3" || "$3" == "-" ]];then
#    bed="/home/pedro.barbosa/mcfonseca/shared/genomes/human/hg38/gencode.v31.annotation.bedops.bed" #fails
    bed="/home/pedro.barbosa/mcfonseca/shared/genomes/human/hg38/gencode.v31.annotation.ucsc.bed12"
else
    bed=$(readlink -f "$3")
fi

if [[ ! -d $(readlink -f "$2") ]]; then
    mkdir $(readlink -f "$2")
fi
OUT=$(readlink -f "$2")

cat > runRSeQC.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=rseqc
#SBATCH --array=0-$(( $JOBS - 1 ))%10
#SBATCH --time=72:00:00
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=mcfonsecalab/rna_quality_control:latest
#SBATCH --output=%j_rseqc.log

readarray -t bams < <(cat $bams)
scratch_out=/home/pedro.barbosa/scratch/rna_seq/\$SLURM_JOB_ID
mkdir \$scratch_out
cd \$scratch_out
srun="srun -N1 -n1"
OUT_BASENAME=\$(basename \${bams[\$SLURM_ARRAY_TASK_ID]})
OUT_BASENAME=\${OUT_BASENAME/Aligned.sortedByCoord.out.bam/}

echo "\$OUT_BASENAME SAMPLE!" 
echo "Inferring type of experiment.."
srun shifter infer_experiment.py -r $bed  -i \${bams[\$SLURM_ARRAY_TASK_ID]} -s 500000 -q 20 > \${OUT_BASENAME}_experiment.txt

echo "Inferring innere distance.."
srun shifter inner_distance.py -r $bed -i \${bams[\$SLURM_ARRAY_TASK_ID]} -s 5 -q 20 -o \${OUT_BASENAME}_inner_distance

echo "Evaluating junction saturation.."
srun shifter junction_saturation.py -r $bed -i \${bams[\$SLURM_ARRAY_TASK_ID]} -o \${OUT_BASENAME}_junction_saturation --percentile-step=10 --min-coverage=2 -q 20

echo "Checking read distribution over genomic features.."
srun shifter read_distribution.py -r $bed -i \${bams[\$SLURM_ARRAY_TASK_ID]} > \${OUT_BASENAME}_read_distribution.txt

echo "Analyzing read duplication levels.."
srun shifter read_duplication.py -i \${bams[\$SLURM_ARRAY_TASK_ID]} -q 20 -o \${OUT_BASENAME}_read_duplication 

echo "Estimating RNA fragment size.."
srun shifter RNA_fragment_size.py -r $bed -i \${bams[\$SLURM_ARRAY_TASK_ID]} > \${OUT_BASENAME}_fragment_size.txt

rm *\.r
mv * $OUT
cd ../ && rm -rf \$SLURM_JOB_ID
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
EOL

sbatch runRSeQC.sbatch

