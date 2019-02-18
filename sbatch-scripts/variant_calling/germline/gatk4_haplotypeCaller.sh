#!/bin/bash

display_usage(){
 printf "Script to run Haplotype caller pipeline in lobo for a bunch of bam files. Be aware that this script already implements the latest innovations regarding
   the multisample variant calling in GVCF mode and ImportGenomicsDB.\n
 Usage:.
    -1st argument must be the file containing the alignment files. A '.bai' index for each bam file must exist in each bam directory.
    -2nd argument must be the name of the final ouptut directory.
    -3rd argument must refer to the type of analysis in hand. Values:[WGS|targeted].
    -4th argument must be the fasta reference for which the reads were aligned. If '-' is set, the existing hg38 version in lobo will be used. Be aware that a fasta index file (.fai) must also be present in the fasta directory.
    -5th argument is optional. If set to false, GNU parallel will be disabled to run the set of samples provided in the 1st argument. Options: [true|false]. Default: true, GNU parallel is used to parallelize the job.
    -6th argument must be provided when the data comes from a targeted experiment. Refers to the target intervals in bed format. Use '-' to skip this argument.
    -7th argument is optional. Refers to the number of nodes,tasks and cpus per task, respectively, to employ on this slurm job in lobo (tasks will be set in parallel,not in the srun command). Default:1,8,5 if GNU parallel (5th argument) is true; 1,1,40 if GNU parallel is false. '-' skips this argument. Setting this argument will override any settings assumed by GNU parallel option.
    -8th argument is optional. Refers to run HaplotypeCaller using SPARK or not. Default:true.
    -9th argument is optional. If set, refers to a known variants file (e.g dbsnp) with IDs.  Its purpose is to annotate our variants with the corresponding reference ID. Must be blocked compressed.\n"
}

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ]; then
    printf "Error. Please set the required parameters of the script\n"
    display_usage
    exit 1
fi

####CHECK BAM INPUT####
if [ ! -f "$1" ]; then
    echo "File "$1" not found!"
    exit 1
fi
while read line
do
    if [ ! -f "$line" ]; then
         printf "Error. $line doesn't exist. Please check more carefully files passed in the 1st argument.\n"
         display_usage
         exit 1
    fi
done < "$1"
BAM_DATA=$(readlink -f "$1")

####SCRATCH WORKDIR####
OUTDIR=$(readlink -f "$2")
WORKDIR="/home/pedro.barbosa/scratch/gatk/haplotypeCaller"
if [ ! -d $WORKDIR ];then
    mkdir $WORKDIR
fi
if [ ! -d $OUTDIR ];then
    mkdir $OUTDIR
fi

##REFERENCE##
if [ "$4" = "-" ]; then
    REF="/mnt/nfs/lobo/IMM-NFS/genomes/hg19/Sequence/WholeGenomeFasta/genome.fa"
elif [ ! -f "$4" ]; then
    printf "Please provide a valid fasta file in th 4th argument.\n"
    display_usage
    exit 1
elif [ ! -f "$(readlink -f ${4}.fai)" ]; then
    printf "Fasta index ${4}.fai not found in the reference directory. Please create one with samtools faidx.\n"
    display_usage
    exit
else
    REF="$(readlink -f "$4")"
fi

re='^[0-9]+$'
##JOB SETTINGS""
PARALLEL="true"
if [[ -z "$7" || "$7" = "-" ]]; then
    if [ -z "$5" ] || [ "$5" = "true" ]; then
        NODES=1
        NTASKS=8
        CPUS=5
    elif [ "$5" = "false" ]; then
        NODES=1
        NTASKS=1
        CPUS=40
        PARALLEL="false"
    else
        printf "Please set a valid value for the 5th argument (gnu parallel)\n"
        display_usage
        exit 1
    fi
else
    IFS=','
    read -r -a array <<< "$7"
    if [ ${#array[@]} = 3 ]; then
        for elem in "${array[@]}"
        do
            if ! [[ "$elem" =~ $re ]]; then 
                printf "Error. Please set INT numbers for the number of nodes, tasks and cpus per task.\n"
                display_usage
                exit 1
            fi
        done
        NODES=${array[0]}
        NTASKS=${array[1]}
        CPUS=${array[2]}
        if [ "$5" = "false" ]; then
            PARALLEL="false"
        elif [ "$5" != "true" ]; then
            printf "Error. Please set a valid value for the 5th argument (gnu parallel flag).\n"
            display_usage
            exit 1
        fi
    else 
        printf "ERROR. 3 fields are required for the 9th argument (nodes,tasks,cpus per task). You set a different number.\n"
        display_usage
        exit 1
    fi
fi

isSpark="true"
if [ -z "$8" ] || [ "$8" == "true" ] || [ "$8" == "-" ]; then
    CMD="gatk ${JAVA_Xmx} HaplotypeCallerSpark"
#    SPARK="-- --spark-runner LOCAL --spark-master local[$CPUS] --num-executors 7 --executor-cores 5 --executor-memory 32500" #single JVM only allows one executor
#    SPARK="-- --spark-runner LOCAL --spark-master local-cluster[7,5,32000]" #simulation of a local pseudo-cluster. failed
    SPARK="-- --spark-runner LOCAL --spark-master local[$CPUS,4]"
#    SPARK="-- --spark-runner SPARK --spark-master"
elif [ "$8" == "false" ]; then
    CMD="gatk ${JAVA_Xmx} HaplotypeCaller"
    isSpark="false"
else
    printf "ERROR. Please set a valid value for the 8th argument. [true|false]\n"
    display_usage
    exit 1
fi

if [ -n "$9" ] && [ -f "$9" ]; then
    printf "Known variants file provided.\n"
    dbsnp=$(readlink -f "$9")
    CMD="$CMD --dbsnp $dbsnp"
elif [ -n "$9" ]; then
    printf "9th argument is not valid file. If you want to use it, please set it properly.\n"
    display_usage
    exit 1
fi

##CMD##
ref_2bit="$(basename $REF)"
#CMD="$CMD --tmp-dir=/home/pedro.barbosa/scratch --reference ${ref_2bit%.*}.2bit -ERC GVCF"
CMD="$CMD --tmp-dir=/home/mcfonseca/scratch --reference $REF -ERC GVCF"
#createOutputVariantIndex true [missing this arg. Needed to add extra step of generating index for genomicsDB import through IndexFeatureFile utility]

###MODE####
analysis=("WGS" "targeted")
BED=$(readlink -f "$6")
if [[ !  " ${analysis[@]} " =~ " ${3} " ]]; then
    printf "Please set a valid value for the 3th argument.\n"
    display_usage
    exit 1
elif [ "$3" != "WGS" ]; then
    if [ -z "$6" ]; then
        printf "When exome or targeted sequencing, you should provide the target regions in the 6th argument as a single bed file.\n"
        display_usage
        exit 1
    elif [[ ${6: -4} == ".bed"  ]]; then
        CMD="$CMD --intervals="$BED" --interval-padding 0"
    else
        printf "Invalid bed input for target regions.\n"
        display_usage
        exit 1
    fi  
else
    CMD="$CMD --intervals='chr1' --interval-padding 0"

# --strict
fi

cat > $WORKDIR/haplotypeCaller_GVCF.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=gatk_hc
#SBATCH --time=72:00:00
#SBATCH --mem=240G
#SBATCH --nodes=$NODES
#SBATCH --ntasks=$NTASKS
#SBATCH --cpus-per-task=$CPUS
#SBATCH --image=broadinstitute/gatk:4.1.0.0
##SBATCH --workdir=$WORKDIR
#SBATCH --output=%j_gatk_hc.log

timestamp() {
    date +"%Y-%m-%d  %T"
}
export -f timestamp
printf "\$(timestamp) -> Job settings\n\$(timestamp) -> isWithParallel=$PARALLEL\n\$(timestamp) -> isSparkjob=$isSpark\n\$(timestamp) -> nodes=$NODES\n\$(timestamp) -> ntasks=$NTASKS\n\$(timestamp) -> cpus=$CPUS\n"
scratch_out=$WORKDIR/\$SLURM_JOB_ID
mkdir \$scratch_out
cd \$scratch_out
srun="srun -N1 -n1 --slurmd-debug 3 shifter"
parallel="parallel --delay 0.2 -j \$SLURM_NTASKS  --env timestamp --joblog parallel.log --halt soon,fail=1 --resume-failed"
echo "\$(timestamp) -> Job started!"
echo "\$(timestamp) -> Converting reference fasta file to 2bit format for HaplotypeCallerSpark"
\$srun --image=mcfonsecalab/htstools_plus:latest faToTwoBit $REF ${ref_2bit%.*}.2bit
echo "\$(timestamp) -> Done"
if [[ -n "$dbsnp"  && ! -f "${dbsnp}.tbi" ]]; then
    echo "\$(timestamp) -> Known variants file was provided, but there is not index. It will be created now.."
    \$srun gatk IndexFeatureFile -F "$dbsnp"
    echo "\$(timestamp) -> Done."
fi

if [ "$PARALLEL" == "true" ]; then
    if [ $isSpark == "true" ]; then
        cat $BAM_DATA | \$parallel  'srun -N1 -n1 --slurmd-debug 3 shifter $CMD -I={} -O={/.}.g.vcf $SPARK && echo -e "{/.}\\t{/.}.g.vcf" >> aux_sampleName.map && shifter gatk IndexFeatureFile -F {/.}.g.vcf'
        #readarray -t bams < $BAM_DATA
        #srun shifter $CMD -I=${bams[$SLURM_ARRAY_TASK_ID]} -O=${bams[$SLURM_ARRAY_TASK_ID]%.*}.g.vcf $SPARK && echo -e "${bams[$SLURM_ARRAY_TASK_ID]%.*}\\t${bams[$SLURM_ARRAY_TASK_ID]%.*}.g.vcf >> aux_sampleName.map && shifter gatk IndexFeatureFile -F ${bams[$SLURM_ARRAY_TASK_ID]%.*}.g.vcf
    else
        cat $BAM_DATA | \$parallel 'srun -N1 -n1 --slurmd-debug 3 shifter $CMD -I={} -O={/.}.g.vcf && echo -e "{/.}\\t{/.}.g.vcf" >> aux_sampleName.map && shifter gatk IndexFeatureFile -F {/.}.g.vcf'
    fi
else
    unique_samples=()
    for j in \$(find $BAM_DATA -exec cat {} \; );do
        printf "\$(timestamp): Processing \$(basename \$j) file.\n"
        i=\$(basename \$j)
        out=\$(echo "\$i" | cut -f1 -d "_")
        if [[ " \${unique_samples[@]} " =~ " \${out} " ]]; then
            printf "\$(timestamp): Warning, duplicate sample ID (\${out}) found after splitting by first '_'. Will use the filename instead without the extension."
            out="\${i%.*}"
        fi
        unique_samples+=(\${out})
        if [ $isSpark == "true" ]; then
            echo $SLURM_JOB_NODELIST
            \$srun $CMD -I=\$j -O=\${out}.g.vcf $SPARK 
        else
            \$srun $CMD -I=\$j -O=\${out}.g.vcf
        fi
        echo -e "\${out}\\t\${out}.g.vcf" >> aux_sampleName.map
        \$srun gatk IndexFeatureFile -F \${out}.g.vcf
        printf "\$(timestamp): Done!\n"
    done
fi
echo "\$(timestamp) -> Haplotype caller in GVCF for all samples has finished."
echo "\$(timestamp) -> Merging gVCF files with CombineGVCFs tools."
CMD_CombineGVCFs="gatk CombineGVCFs -R $REF -O combined.g.vcf.gz"
#echo "\$(timestamp) -> Creating GenomicsDB datastore to combine multiple gvcf data."
##TRY ALSO FOR ONE SINGLE INTERVALS FILE###
#CMD_GenomicsDB="gatk GenomicsDBImport --sample-name-map aux_sampleName.map --batch-size 5 -R $REF --genomicsdb-workspace-path workspace"
while read line;do
    gvcf_file=\$(echo "\$line" | cut -f2)
    CMD_CombineGVCFs="\$CMD_CombineGVCFs -V \$gvcf_file"
done < aux_sampleName.map
if [ "$3" == "WGS" ]; then
    \$srun \$CMD_CombineGVCFs
#    for ((i=1;i<=22;i++)); do echo -e "\$(timestamp) Consolidating gVCF for chromosome \${i}"; \$srun \$CMD_GenomicsDB --L chr\${i}; done
#    echo -e "\$(timestamp) Consolidating gVCF for sexual chromosomes and mithocondria."
#    \$srun $CMD_GenomicsDB -L chrM
#    \$srun $CMD_GenomicsDB -L chrY
#    \$srun $CMD_GenomicsDB -L chrX
else
    \$srun \$CMD_CombineGVCFs
#    while read line;
#    do
#        echo -e "\$(timestamp) Consolidating gVCF for \$(echo "\$line" | cut -f1,2,3) interval."
#        \$srun \$CMD_GenomicsDB --L \$(echo "\$line" | cut -f1):\$((\$(echo "\$line" | cut -f2) + 1))-\$(echo "\$line" | cut -f3)
#    done < $BED
fi
echo "\$(timestamp) -> Done! Genotyping gVCFs stored in the GenomicsDB database."
\$srun gatk GenotypeGVCFs --variant combined.g.vcf.gz -R $REF --output joint.vcf.gz
#\$srun gatk GenotypeGVCFs --variant gendb://workspace -R $REF --output joint.vcf.gz
mv * ../\${SLURM_JOB_ID}_gatk_hc.log $OUTDIR
cd .. && rm -rf \$SLURM_JOB_ID

echo -e "\$(timestamp) -> Finished job."
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
echo -e "\$(timestamp) -> All done!"
EOL
sbatch $WORKDIR/haplotypeCaller_GVCF.sbatch
sleep 3 

cd $WORKDIR
mv haplotypeCaller_GVCF.sbatch $(ls -td -- */ | head -n 1)
