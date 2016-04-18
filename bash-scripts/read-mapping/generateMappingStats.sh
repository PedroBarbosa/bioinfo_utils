#!/bin/bash
usage="$(basename "$0") input_file mapper [-h] [-b] [-s]. Script to generate different mapping statistics.

where:
    input_file	Text file where each line represents the path for the SAM/BAM files to be processed.
    mapper	Aligner used to generate the alignments. Available options: Tophat2, Hisat2, STAR, bowtie2.
   
    Optional arguments:
    -h  show this help text
    -b  flag to indicate that the bam/sam files are sorted by read names. Might be useful to gain speed when calculating unique mapped reads.
    -s  flag to indicate that the reads are single-end. (default paired-end)    "

seed=42
while getopts ':hs:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    s) seed=$OPTARG
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND - 1))

if [ -z $1 ] || [ -z $2 ]; then
    printf "ERROR: Please set all the required arguments for the script.\n\n"
    echo "$usage" >&2
    exit 1
fi


containsElement () {
  local e
  for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
  return 1
}

MAPPERS_AVAILABLE=("Tophat2" "Hisat2" "STAR","bowtie2")
containsElement "$2" "${MAPPERS_AVAILABLE[@]}"
if [ $(echo $?) = 1 ]; then
    printf "ERROR: Unrecognized mapper. Please set a valid aligner name.\n\n"
    echo "$usage" >&2
    exit 1	
else
    MAPPER=$2
fi

#ARG1=${@:$OPTIND:1}
#ARG2=${@:$OPTIND+1:1}


if [ -f $PWD/mappingStats.txt ]; then
	rm $PWD/mappingStats.txt

fi

command -v samtools >/dev/null 2>&1 || { echo >&2 "Samtools is required to be on the system path, but it seem that it's not installed.Aborting."; exit 1; }
while read line
do
    filename=$line
    printf "Generating mapping stats for $filename file..\n" 2>&1 | tee -a $PWD/mappingStats.txt
    printf '\n%s\n' '####FINAL STATS####' >> $PWD/mappingStats.txt
    printf '%s\t%s\n' 'Number of alignments:' $(samtools view -c $filename) >> $PWD/mappingStats.txt
    printf '%s\t%s\n' 'Number of unmapped reads:' $(samtools view -cf4 $filename) >> $PWD/mappingStats.txt
    printf '%s\t%s\n' 'Number of primary and linear aligments:' $(samtools view -cF2308 $filename) >> $PWD/mappingStats.txt
    printf '%s\t%s\n' 'Number of primary and linear aligments with mapping quality > 10:' $(samtools view -cF2308q10 $filename) >> $PWD/mappingStats.txt
    printf '%s\t%s\n' 'Number of secondary alignments:' $(samtools view -cf256 $filename) >> $PWD/mappingStats.txt
    printf '%s\t%s\n' 'Number of chimeric alignments:' $(samtools view -cf2048 $filename) >> $PWD/mappingStats.txt
    printf '%s\t%s\n' 'Number of reads mapped as proper pair:' $(samtools view -cf2 $filename) >> $PWD/mappingStats.txt
    printf '\t%s\t%s\n' 'Number of proper pairs mapped as FR:' $(samtools view -cf99 $filename) >> $PWD/mappingStats.txt
    printf '\t%s\t%s\n' 'Number of proper pairs mapped as RF:' $(samtools view -cf83 $filename) >> $PWD/mappingStats.txt
    if [ $MAPPER = "STAR" ];then
	printf '%s\t%s\n' 'Number of alignments in STAR with the NH:i:1 tag (supposely represent unique mappers):' $(samtools view $filename | grep -wc "NH:i:1") >> $PWD/mappingStats.txt
        printf '%s\t%s\n' 'Number of alignments in STAR with the mapping quality of 255 (supposely represent unique mappers, if you did not change this parameter when running star):' $(samtools view -cq255 $filename) >> $PWD/mappingStats.txt
	printf '%s\t%s\n' 'Number of unique pairs using bitflag, tags and sort/uniq filtering :' $(samtools view -q255 -F 2308 $filename | grep -w "NH:i:1" | cut -f 1 | sort | uniq -c | sed -e 's/^[ \t]*//' | grep "^2" | cut -d ' ' -f 2 | wc -l) >> $PWD/mappingStats.txt
    elif [ $MAPPER = "bowtie2" ];then
        printf '%s\t%s\n' 'Number of alignments in bowtie2 with the mapping quality of 255 (supposely represent unique mappers):' $(samtools view -cq255 $filename) >> $PWD/mappingStats.txt

    fi
    printf "DONE!!\n\n\n" 2>&1 | tee -a $PWD/mappingStats.txt
done <$1
