#!/bin/bash
#SBATCH --job-name=gnomAD_removeFileds
#SBATCH --time=72:00:00
#SBATCH --array=0-23%5
#SBATCH --mem=240G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
##SBATCH --workdir=/home/pedro.barbosa/scratch
#SBATCH --output=%j_gnomAD_removeFileds.log
#SBATCH --image=mcfonsecalab/variantutils:0.4

infile=$(readlink -f $1)
chrFile=$(readlink -f $2)

readarray -t chrs < $chrFile
cd /home/pedro.barbosa/scratch

srun shifter bcftools annotate -r ${chrs[$SLURM_ARRAY_TASK_ID]} -x "^INFO/AC_raw,INFO/AN_raw,INFO/AF_raw,INFO/nhomalt_raw,INFO/AC,INFO/AN,INFO/AF,INFO/nhomalt,INFO/AC_afr,INFO/AN_afr,INFO/AF_afr,INFO/nhomalt_afr,INFO/AC_eas,INFO/AN_eas,INFO/AF_eas,INFO/nhomalt_eas,INFO/AC_amr,INFO/AN_amr,INFO/AF_amr,INFO/nhomalt_amr,INFO/AC_nfe,INFO/AN_nfe,INFO/AF_nfe,INFO/nhomalt_nfe,INFO/AC_nfe_nwe,INFO/AN_nfe_nwe,INFO/AF_nfe_nwe,INFO/nhomalt_nfe_nwe,INFO/AC_nfe_seu,INFO/AN_nfe_seu,INFO/AF_nfe_seu,INFO/nhomalt_nfe_seu" -Oz -o chr_${chrs[$SLURM_ARRAY_TASK_ID]}_gnomAD.vcf.gz $infile

mv chr_${chrs[$SLURM_ARRAY_TASK_ID]}_gnomAD.vcf.gz $(dirname $infile)

