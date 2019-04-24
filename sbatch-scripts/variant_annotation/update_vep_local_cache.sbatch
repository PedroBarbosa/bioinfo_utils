#!/bin/bash
display_usage(){
    printf "
This is the args allowed, all of them are optional.
1st argument is the species name. Values: [homo_sapiens,mus_musculus,-]. Default: homo_sapiens
2nd argument is the assembly name. Values: [GRCh38, GRCh37, -]. Default: GRCh38.
3rd argument is the cache type to install: [refseq,merged,-]. Default: ensembl.
4rd argument is the cache version. Default: latest.\n"
}

display_usage
if [[ -z "$1" || "$1" == "-" ]]; then
    SPECIES="homo_sapiens"  
elif [[ "$1" == "mus_musculus" ]]; then
    SPECIES="mus_musculus"
else
    printf "Please select a valid value for the 1st argument.\n"
    display_usage
    exit 1
fi      


if [[ -z "$2" || "$2" == "-" ]]; then
    if [[ "$SPECIES" == "homo_sapiens" ]]; then
        ASSEMBLY="GRCh38"
    else
        ASSEMBLY="GRCm38"
    fi  
elif [[ "$2" == "GRCh37" ]]; then
    ASSEMBLY="GRCh37"
else
    printf "Please select a valid value for the 2nd argument.\n"
    display_usage
    exit 1
fi      


if [[ -z "$3" || "$3" == "-" ]]; then
    CACHE=${SPECIES}  
elif [[ "$3" == "refseq" ]]; then
    CACHE="${SPECIES}_refseq"
elif [[ "$3" == "merged" ]]; then
    CACHE="${SPECIES}_merged"
else
    printf "Please select a valid value for the 3rd argument.\n"
    display_usage
    exit 1
fi     

if [[ ! -z "$4" ]]; then
    CACHE_VERSION="$4"
fi  

CACHE_DIR="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/cache"
cat > $PWD/update_vep_local_cache.sbatch <<EOL
#!/bin/bash
#SBATCH --job-name=updateCache
#SBATCH --time=72:00:00
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
##SBATCH --workdir=/home/pedro.barbosa/scratch/vep
#SBATCH --output=%j_vep_cache.log

mkdir /home/pedro.barbosa/scratch/vep/$SLURM_JOB_ID && cd /home/pedro.barbosa/scratch/vep/$SLURM_JOB_ID

if [[ -z "$CACHE_VERSION" ]]; then
    srun shifter -V=$CACHE_DIR:/media --image=ensemblorg/ensembl-vep:latest INSTALL.pl -c /media -a cfp -s $SPECIES -y $ASSEMBLY -g all
else
    srun shifter -V=$CACHE_DIR:/media --image=ensemblorg/ensembl-vep:latest INSTALL.pl -c /media -a cfp -s $SPECIES -y $ASSEMBLY -g all --CACHE_VERSION $CACHE_VERSION
fi
EOL

sbatch update_vep_local_cache.sbatch
