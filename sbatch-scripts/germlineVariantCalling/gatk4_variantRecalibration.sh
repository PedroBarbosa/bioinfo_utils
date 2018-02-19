#!/usr/bin/env bash
display_usage(){
 printf "Script to run GATK4 Variant Recalibration pipeline in lobo.\n
 Usage:.
    -1st argument must be the input VCF file to recalibrate.
    -2nd argument must be the name of the final ouptut directory.
    -3rd argument must refer to the type of analysis in hand. Values:[WGS|targeted]. Useful to apply different values in some arguments.
    -4th argument must be the fasta reference for which the reads were aligned. If '-' is set, the existing hg38 version in lobo will be used.
     Be aware that a fasta index file (.fai) must also be present in the fasta directory.
    -5th argument must be a file listing the resources to employ for in the recalibration. Must be a 2 column tab separated file where first
    column must refer to the input file and second to the resource name. Available resources names: [hapmap,omni,1000G,dbsnp,mills]. At least
     one resource is mandatory.
    -6th argument is optional. Refers to the names of the annotations that should used for calculations. Use '-' to skip the argument. Multiple
    annotations must be splitted by a comma. Default annotations to be used:[
    -5th argument is optional. If set to false, GNU parallel will be disabled to run the set of samples provided in the 1st argument. Options: [true|false]. Default: true, GNU parallel is used to parallelize the job.
    -6th argument must be provided when the data comes from a targeted experiment. Refers to the target intervals in bed format. Use '-' to skip this argument.
    -7th argument is optional. If set, refers to a known variants file (e.g dbsnp) with IDs.  Its purpose is to annotate our variants with the corresponding reference ID.\n"
}

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ]; then
    printf "Error. Please set the required parameters of the script\n"
    display_usage
    exit 1
fi

