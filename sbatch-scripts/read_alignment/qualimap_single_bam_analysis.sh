#!/bin/bash
display_usage(){
 printf "Script to automatically run qualimap quality control analysis for a set of bam files.\n
 Usage:
    -1st argument must the bam files to process.
    -2nd argument must the output directory.
    -3th argument must the annotation in GTF format. Default (hg38 gencode v33 primary assembly). Use '-' to skip this argument.
    -4th argument must the strandness protocol of the data. Default (strand-specific-reverse). Possible values: (strand-specific-forward, non-strand-specific). Use '-' to skip this argument.
    -5th argument must be the species under analysis. Default: HUMAN. Possible values: (MOUSE).Use '-' to skip this argument.\n"
}

if [ -z "$1" ] || [ -z "$2" ]; then
    printf "Error. Please set the required parameters\n"
    display_usage
    exit 1
fi

bams=$(readlink -f "$1")
JOBS=$(cat $bams | wc -l)

if [[ ! -d $(readlink -f "$2") ]]; then
    mkdir $(readlink -f "$2")
fi
OUT=$(readlink -f "$2")

if [[ -z "$3" || "$3" == "-" ]]; then
    GTF="/home/pedro.barbosa/mcfonseca/shared/genomes/human/hg38/gencode.v33.primary_assembly.annotation.gtf"
else
    GTF=$(readling -f "$3")
fi

if [[ -z "$4" || "$4" == "-" || "$4" == "strand-specific-reverse" ]]; then
    STRANDNESS="strand-specific-reverse"
elif [[ "$4" == "strand-specific-forward" || "$4" == "non-strand-specific" ]]; then
    STRANDNESS="$4"
else
    printf "Please set a valid argument for the 4th argument.\n"
    display_usage
    exit 1
fi

if [[ -z "$5" || "$5" == "-" || "$5" == "HUMAN" ]];then
    SPECIES="HUMAN"
elif [[ "$5" == "MOUSE" ]]; then
    SPECIES="$5"
else
    printf "Please set a valid argument for the 5th argument.\n"
    display_usage
    exit 1
fi  

cat > qualimap.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=qualimap
#SBATCH --array=0-$(( $JOBS - 1 ))%5
#SBATCH --time=72:00:00
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --image=pegi3s/qualimap:latest
#SBATCH --output=%j_qualimap.log

readarray -t bams < <(cat $bams)
scratch_out=/home/pedro.barbosa/scratch/rna_seq/\$SLURM_JOB_ID
mkdir \$scratch_out
cd \$scratch_out
srun="srun -N1 -n1"
OUT_BASENAME=\$(basename \${bams[\$SLURM_ARRAY_TASK_ID]})
OUT_BASENAME=\${OUT_BASENAME/bam/}

echo "\$OUT_BASENAME SAMPLE!"
####### General BAM QC #######
#mkdir bamqc 
#cd bamqc
#srun shifter qualimap bamqc --java-mem-size=8G -bam \${bams[\$SLURM_ARRAY_TASK_ID]} -outdir . --genome-gc-distr $SPECIES \
#-gff $GTF --collect-overlap-pairs \
#-nt \$SLURM_CPUS_PER_TASK \
#--outside-stats \
#--sequencing-protocol $STRANDNESS
#cd ../

######RNA QC#######
mkdir rnaqc
cd rnaqc
srun shifter qualimap rnaseq --java-mem-size=8G -bam \${bams[\$SLURM_ARRAY_TASK_ID]} -gtf $GTF --sequencing-protocol $STRANDNESS -outdir \$PWD
cd ../../

mv \$SLURM_JOB_ID \$OUT_BASENAME 
mv \$OUT_BASENAME $OUT
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
EOL

sbatch qualimap.sbatch

