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
cache_dir="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/cache"

srun shifter -V=$cache_dir:/media --image=ensemblorg/ensembl-vep:latest INSTALL.pl -c /media -a cfp -s homo_sapiens, mus_musculus -y GRCh38 -g AncestralAllele,Blosum62,CADD,CSN,Carol,Condel,Conservation,Downstream,Draw,ExAC,ExACpLI,FATHMM_MKL,G2P,GO,GeneSplicer,Gwava,LD,LOVD,LoF,LoFtool,LocalID,MPC,MTR,MaxEntScan,NearestGene,Phenotypes,ProteinSeqs,REVEL,SameCodon,SpliceRegion,TSSDistance,dbNSFP,dbscSNV,miRNA

#srun shifter -V=$cache_dir:/media --image=ensemblorg/ensembl-vep:latest INSTALL.pl -c /media -a cfp -s homo_sapiens -y GRCh37 -g AncestralAllele,Blosum62,CADD,CSN,Carol,Condel,Conservation,Downstream,Draw,ExAC,ExACpLI,FATHMM_MKL,G2P,GO,GeneSplicer,Gwava,LD,LOVD,LoF,LoFtool,LocalID,MPC,MTR,MaxEntScan,NearestGene,Phenotypes,ProteinSeqs,REVEL,SameCodon,SpliceRegion,TSSDistance,dbNSFP,dbscSNV,miRNA

#srun shifter -V=$cache_dir:/media --image=ensemblorg/ensembl-vep:latest INSTALL.pl -c /media -a cf -s mus_musculus
