#!/bin/bash


display_usage(){
	printf "First argument must be the file of contigs to scaffold.
Second argument must a file that lists the path where all the sorted and indexed bam files are located. Files must be sorted by increasing insert size value, it is a besst requirement.
Third argument must be the orientation of the libraries aligned in the bam file. Available option: [fr|rf].
Fourth argument must be the path to the output.
Fifth argument is optional. It refers to the option of providing the insert size and stdv for the mate pair libraries. Available options: [true|false]. Default: false.\n"
}


if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ]; then
	printf "Too few arguments.\n"
	display_usage
	exit 1
fi

CONTIG_FILE="-c $1"

#declare -A list_sorted
#while read line; do
#  list_sorted+=('$line')
#done < $2
declare -a list_sorted
readarray -t list_sorted < $2 # Include newline.


BAM_FILES="-f ${list_sorted[@]}"
ORIENTATION=""
OUTPUT="-o $4"
INSERT_SIZES_MP="-m 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 2100 5000 5000 5000 5000 5000 5000"
INSERT_STDV_MP="-s 250 250 250 250 250 250 250 250 250 250 250 250 600 600 600 600 600 600"
#for loop to givee orientation for each library
for i in "${list_sorted[@]}"
do
#	BAM_FILES="$BAM_FILES $i"
	ORIENTATION="$ORIENTATION$3 "
done

ORIENTATION="--orientation $ORIENTATION"

if [ -z "$5" ] || [ "$5" = "false" ]; then
	command="runBESST -plots -q --min_mapq 10 $CONTIG_FILE $BAM_FILES $OUTPUT $ORIENTATION"
elif [ "$5" = "yes" ]; then
	command="runBESST -plots -q --min_mapq 10 $CONTIG_FILE $BAM_FILES $OUTPUT $ORIENTATION $INSERT_SIZES_MP $INSERT_STDV_MP"
else
	printf "Not a valid option for the mate pair manual insert size flag [5th argument]"
	display_usage
	exit 1
fi

##Pass command to script and give permissions to run
file="runBESST.sh"
echo -e "#!/bin/bash\n$command 2> error.log" > ./$file | chmod 755 $file
