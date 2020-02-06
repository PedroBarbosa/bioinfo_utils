#!/bin/bash
#This is a script that runs the insert size estimation script from SSPACE for a file containing the path for a list of different mate pair  ordered pairs
ispair1=true
ispair2=false
matepairFlag=false
num_samples=0
contigs_file=/mnt/msa/assembly/ray-assembly/highstringency/PE/k51/output/Contigs.fasta
cd /opt/tools/SSPACE-STANDARD-3.0_linux-x86_64/tools
while read line
do
        filename=$line
        if [[ "$filename" == *"MP"* ]] ; then
                matepairFlag=true
        fi
        if [[ "$matepairFlag" = "true" && "$ispair1" = "true" && "$ispair2" = "false" ]]; then
                file1=$filename
                ispair1=false
                ispair2=true
        elif [[ "$matepairFlag" = "true" && "$ispair1" = "false" && "$ispair2" = "true" ]]; then
                file2=$filename
                ispair1=true
                ispair2=false
                echo "Processing..."
                let "num_samples +=1"
                perl estimate_insert_size.pl -c $contigs_file -1 $file1 -2 $file2 -m 6000 -p testinsertMP > MP${num_samples}.txt
                echo $num_samples processed..
        fi
done < $1
