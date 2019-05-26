#!/usr/bin/env bash
display_usage() {
echo 'Script to run Hisat2 for multiple fastq files.

-1st argument must be the file listing RNA-seq pairs consecutively. One file per line.
-2nd argument must be the directory of the reference indexed database. Default: human hg38
-3rd argument must be the output directory.
-4th argument is optional. It is the identifier to extract the sample pair names from fastq files. Default: "_1.fastq".'
}

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] ; then
        printf "ERROR:Please provide at least the 3 first arguments required for the script.\n\n"
        display_usage 
        exit 1
fi

readarray FASTQ < $(readlink -f "$1")
JOBS=$(( ${#FASTQ[@]} / 2 ))
if [[ $2 == "-" ]]; then
    INDEX="/home/pedro.barbosa/mcfonseca/shared/genomes/human/hg38/hisat2/hg38_hisat2"
else
    INDEX=$(readlink -f "$2")
fi


if [[ ! -d $(readlink -f "$3") ]]; then
    mkdir $(readlink -f "$3")
fi
OUT=$(readlink -f "$3")

if [[ -z "$4" ]]; then
    DEL="_1.fastq"
    p2="_2.fastq"
else
    DEL="$4"
    p2=${DEL/1/2}
fi

cat > hisat2.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=hisat2
#SBATCH --array=0-$(($JOBS - 1))%5
#SBATCH --time=72:00:00
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --image=mcfonsecalab/hisat2pipeline:latest
#SBATCH --output=%j_hisat2.log

readarray -t pair1 < <(cat \$(readlink -f "$1") | grep -e "R1.fastq" -e "R1.fq" -e "_1.fastq" -e "_1.fq")
scratch_out=/home/pedro.barbosa/scratch/rna_seq/hisat2/\$SLURM_JOB_ID
mkdir \$scratch_out
cd \$scratch_out
srun="srun -N1 -n1"

dir=\$(dirname \${pair1[\$SLURM_ARRAY_TASK_ID]})
pair2=\$(basename \${pair1[\$SLURM_ARRAY_TASK_ID]/$DEL/$p2})
fullpathpair2="\$dir/\$pair2"

outfile=\${fullpathpair2##*/}
outfile=\${outfile%.*}
outfile=\${outfile/$p2/}.bam

srun shifter --image=docker:mcfonsecalab/hisat2pipeline:latest hisat2 --rna-strandness RF --fr --no-discordant --no-mixed --summary-file \${outfile/.bam/_summary.txt} -p \$SLURM_CPUS_PER_TASK  -x $INDEX -1 \${pair1[\$SLURM_ARRAY_TASK_ID]} -2 \$fullpathpair2 | shifter --image=docker:mcfonsecalab/htstools_plus:latest samtools view -Shu - | shifter --image=docker:mcfonsecalab/htstools_plus:latest samtools sort -O bam -@ \$CPUS_PER_TASK - > \$outfile 
srun shifter --image=docker:mcfonsecalab/htstools_plus:latest samtools index \$outfile

echo "Extracting uniques.."
srun shifter --image=docker:mcfonsecalab/htstools_plus:latest samtools view -h -f2 \$outfile | awk 'substr(\$1, 0, 1)=="@" || \$0 !~ /ZS:/' | shifter --image=docker:mcfonsecalab/htstools_plus:latest samtools view -hb > \${outfile/.bam/_uniq.bam}
#srun shifter --image=docker:mcfonsecalab/htstools_plus:latest samtools index \${outfile/.bam/uniq.bam}
mv * $OUT

#Single-end Human
#$parallel  '$srun  shifter --image=docker:rluis/hisat2pipeline:1.0  hisat2 --rna-strandness RF --fr --no-discordant  --no-mixed --summary-file summary_file.txt -p 40 -x /home/rluis/Rui-testing/Genome/hg38_HISAT2_INDEX/Ensembl/ensemblFasta/hg38_EnsemblGRCh38_Ensembl90  -U {1} |  shifter --image=docker:mcfonsecalab/htstools_plus:latest  samtools view -@ 40 -hf 3  -o Alignment.sam - ' ::: $full_name1

#Filter Unique Reads
#grep  -e "^@*" -we "NH:i:1"  Alignment.sam > Alignment_Unique.sam

cd ../ && rm -rf \$SLURM_JOB_ID
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
rm -rf `ls -la /tmp/ | grep 'pedro.barbosa' | awk ' { print $9 } '`
EOL

sbatch  hisat2.sbatch
