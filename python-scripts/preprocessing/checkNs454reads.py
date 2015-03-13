__author__ = 'pedro'

import argparse
import os
import csv
import collections

def process_fasta_file(fastafile):
    map={}
    number_of_Ns = 0
    number_of_reads = 0
    first_read = True
    read = ""
    reads_Ns_at_end = 0
    with open(fastafile[0]) as file:
        for line in file:
            line=line.rstrip().upper();
            if line.startswith('>') and first_read:
                number_of_reads += 1
                first_read = False
            elif line.startswith('>'):

                read_length = len(read)
                read_NoN = len(read.rstrip("N"))
                if read_length != read_NoN:
                    reads_Ns_at_end += 1


                if number_of_Ns in map:
                    map[number_of_Ns] += 1
                else:
                    map[number_of_Ns] = 1


                number_of_reads += 1
                number_of_Ns = 0
                read= ""
            else:
                read += line
                number_of_Ns += line.count("N")


        #process last read
        if(number_of_Ns) in map:
            map[number_of_Ns] += 1
        else:
            map[number_of_Ns] = 1

    print("Number of reads:\t" ,number_of_reads)
    print("Reads with Ns at the end",reads_Ns_at_end)
    file.close()
    return map, number_of_reads




def write2file(dict, numb_reads, outputfile):

    if os.path.exists(outputfile):
        os.remove(outputfile)
    with open(outputfile, "w") as csvfile:
        writer = csv.writer(csvfile,dialect=csv.excel_tab)
        writer.writerow(('Total number of reads:', numb_reads))
        writer.writerow('')
        writer.writerow(('Number N\'s in reads', 'Number of reads with given number of N\'s', 'Percentage (%)'))


        od = collections.OrderedDict(sorted(dict.items()))
        for k, v in od.iteritems():
            percentage = round((v / float(numb_reads) * 100), 3)
            writer.writerow((k,v,percentage))

        csvfile.close()


parser = argparse.ArgumentParser(description='Script to analyze the pattern of Ns in a fasta file of reads.')
parser.add_argument(dest='fasta_file', metavar='fasta_file', nargs=1,
            help='File of reads to be processed.')
parser.add_argument(dest='output_file', metavar='output_file', nargs=1, help='Output file')
args = parser.parse_args()


if __name__ == "__main__":
    print("Counting Ns in reads..")
    dict, numb_reads = process_fasta_file(args.fasta_file)
    print("Writing the output to " + args.output_file[0] + "..")
    write2file(dict, numb_reads, args.output_file[0])
    print("Done.")
