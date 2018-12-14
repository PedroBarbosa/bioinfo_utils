#!/bin/bash
display_usage(){
 printf "Script to automatically run bwa mem in lobo for multiple Illumina paired-end files.\n
 Usage:.
    -1st argument must be the file containing the sequences to align. Each pair needs come consecutively.
    -2nd argument must be the name of the ouptut directory. If '-' is given, output will be written in the parent directory of the fastq files.
    -3rd argument is optional. If set, it refers to the indexed database. By default, human hg19 version is set. If you want to use another reference, or another human version, please set here the path accordingly. If you set '-' here, this parameter will be ignored, and defaults will be employed.
    -4th argument is optional. If set to false, GNU parallel will be disabled to run the set of samples provided in the 1st argument. Options: [true|false]. Default: true, GNU parallel is used to parallelize the job.
    -5th argument is optional. Read groups are added anyway, but if set, this argument must be the string to add to the LB, PL and PU fields, respectively, comma separated. Ex: LB:lib,PL:instrument,PU:platformUnit. The SM and ID fields will be automatically added based on the sample basename (up to the first '_' found), regardless this argument is set or not.Defaults:'LB:pe,PL:hiseq,PU:4000'\n"
}

if [ -z "$1" ] || [ -z "$2" ]; then
    printf "Error. Please set the required parameters of the script\n"
    display_usage
    exit 1
fi

####CHECK FASTQ INPUT####
if [ ! -e "$1" ]; then
    printf "File provided in the 1st argument does not exist. Please set a valid one.\n"
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
FASTQ_DATA=$(readlink -f "$1")

####SCRATCH WORKDIR####
WORKDIR="/home/pedro.barbosa/scratch/bwa"
if [ ! -d $WORKDIR ];then
    mkdir $WORKDIR
fi

####OUTDIR####
if [ "$2" = "-" ]; then
    OUTDIR="{//}"
else
    OUTDIR=$(readlink -f "$2")
    if [ ! -d $OUTDIR ];then
        mkdir $OUTDIR
    fi
fi

####INNDEX####
if [ -z "$3" ] || [ "$3" = "-" ]; then
    INDEX="/mnt/nfs/lobo/IMM-NFS/genomes/hg19/Sequence/BWAIndex/genome.fa"
elif [ -d $(dirname "$3") ]; then
    INDEX="$3"
else
    printf "Error. Please set a valid value for the genome index.\n"
    display_usage
    exit 1    
fi

###PARALLEL####
if [ -z "$4" ] || [ "$4" = "true" ]; then
    NODES=1 #2 #1
    NTASKS=4 #10 #4
    CPUS=10 #8 #10
    PARALLEL=true
elif [ "$4" = "false" ]; then
    NODES=1
    NTASKS=1
    CPUS=40
    PARALLEL=false
else
    printf "Please set a valid value for the 4th argument\n"
    display_usage
    exit 1
fi

##READGROUPS####
if [ -z "$5" ]; then
    READGROUP="LB:pe\tPL:hiseq\tPU:4000"
elif [[ "$5" == *"SM"* || "$5" == *"ID"* ]]; then
    printf "SM and ID read group fields are automatically added. Please just set LB,PL and PU fields.\n"
    display_usage
    exit 1
elif [[ "$5" == *"LB"* && "$5" == *"PL"* && "$5" == *"PU"* && "$5" == *","* ]]; then
    READGROUP=$(echo "${5//,/\\t}")
else
    printf "Please set a valid read group string. Include LB,PL and PU fields and don't forget to separate such fields with a comma.\n"
    display_usage
    exit 1
fi
echo $READGROUP

cat > $WORKDIR/runBwa.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=m2sb
#SBATCH --time=72:00:00
#SBATCH --mem=120G
#SBATCH --nodes=$NODES
#SBATCH --ntasks=$NTASKS
#SBATCH --cpus-per-task=$CPUS
#SBATCH --image=docker:ummidock/bwa
#SBATCH --workdir=$WORKDIR
#SBATCH --output=$WORKDIR/%j_bwa.log

SCRATCH_OUTDIR="$WORKDIR/\$SLURM_JOB_ID"
mkdir \$SCRATCH_OUTDIR
cd \$SCRATCH_OUTDIR
echo "`date`: Analysis started."
if [ "$PARALLEL" = "true" ]; then
	srun="srun --exclusive -N1 -n1"
	parallel="parallel --tmpdir \$SCRATCH_OUTDIR --halt soon,fail=1 --delay 0.2 -j $NTASKS --joblog parallel.log --resume-failed"	
        #Test instead of cat initial fofn file: [parallel --pipepart -a $FASTQ_DATA] and [cat $FASTQ_DATA | parallel --pipe]
        cat "$FASTQ_DATA" | grep -v "R2" | \$parallel '\$srun shifter bwa mem -R "@RG\tID:{=s{.*/}{};s/\_[^_]+$//;=}_id\tSM:{=s{.*/}{};s/\_[^_]+$//;=}\t$READGROUP" -t $CPUS $INDEX {} {=s/R1/R2/=} | shifter --image=ummidock/bowtie2_samtools:latest samtools sort -O bam -@ $CPUS -o {=s{.*/}{};s/\_[^_]+$//;=}_sort.bam -; shifter --image=ummidock/bowtie2_samtools:latest samtools index {=s{.*/}{};s/\_[^_]+$//;=}_sort.bam; mv {=s{.*/}{};s/\_[^_]+$//;=}_sort.bam* $OUTDIR'
else
	for j in \$(find $FASTQ_DATA -exec cat {} \; | grep -v "R2" );do
	    printf "\$(date): Processing \$(basename \$j) file.\n"
	    i=\$(basename \$j)
	    R2=\$(echo \$j | sed 's/R1/R2/')
	    outbam=\$(echo \$i | sed 's/_R1.\+//')
	    sampleID=`echo \$i | cut -f1 -d _`
	    readGroup="@RG\tID:${sampleID}_id\tSM:$sampleID\t$READGROUP"
	    srun shifter bwa mem -R "\$readGroup" -t $CPUS $INDEX \$j \$R2 | shifter --image=ummidock/bowtie2_samtools:latest samtools sort -O bam -@$CPUS - > \${outbam}_sort.bam
	    srun shifter --image=ummidock/bowtie2_samtools:latest samtools index \${outbam}_sort.bam
            if [ $OUTDIR = "{//}" ]; then
                mv \${outbam}_sort.bam* \$(dirname \$j)
            else
                mv \${outbam}_sort.bam* $OUTDIR
            fi
	done
fi 
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID 
if [ $OUTDIR != "{//}" ]; then
    mv \$SCRATCH_OUTDIR/* $OUTDIR
else
    echo "Since you set the output directory to be different for each sample, the run logs weren't moved anywhere. You can find them in the $WORKDIR folder."
fi
echo "`date`: All done."
EOL
sbatch $WORKDIR/runBwa.sbatch
sleep 1
cd $WORKDIR
mv runBwa.sbatch $(ls -td -- */ | head -n 1) #moves sbatch to the last directory created in bwa dir (we expect to be the slurm job ID folder)
