# Arguments:
#   $1: fofn
#
# Return code:
#   0 if at least one line was read.
#   1 if input is empty.

# Description:
#   Reads a fofn by a defining a movie as the set of 3 bax files and runs bax2bam on each chunk
function readlines () {
    local N="$1"
    local line
    local rc="1"

    # Read at most N lines
    for i in $(seq 1 3)
    do
        # Try reading a single line
        read line
        if [ $? -eq 0 ]
        then
            # Output line
            echo $line
            rc="0"
        else
            break
        fi
    done

    # Return 1 if no lines where read
    return $rc
}

while chunk=$(readlines )
do
    bax2bam $chunk --subread --pulsefeatures=DeletionQV,DeletionTag,InsertionQV,IPD,MergeQV,SubstitutionQV,PulseWidth,SubstitutionTag &> log_bax2bam.txt
done < $1
