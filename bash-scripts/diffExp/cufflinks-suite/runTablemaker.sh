#!/bin/bash
while read line
do
SAMPLE=$(basename $line)
printf "Runnning $SAMPLE sample..\n"
/opt/tools/tablemaker-2.1.1.Linux_x86_64/tablemaker -p 50 --library-type fr-firststrand -m 150 -o tableMaker-${SAMPLE} -W -G mergedAssembly/merged.gtf $line
printf "Finished running tablemaker for $SAMPLE sample."
done < $1
printf "DONE\n"
