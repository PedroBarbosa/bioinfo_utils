__author__ = 'pedro'

import argparse
import os
import csv
import collections

def process_qual_file(qualfile):
    map={}
    number_of_reads = 0
    number_of_nucleotides = 0
    with open(qualfile[0]) as file:
        for line in file:
            line=line.rstrip()
            if line.startswith('>'):
                number_of_reads += 1
            else:
                qualities = line.split()
                for val in qualities:
                    number_of_nucleotides += 1
                    if val in map:
                        map[val] += 1
                    else:
                        map[val] = 1

    file.close()
    return map, number_of_nucleotides, number_of_reads


def write2file(dict, numb_nuc, numb_reads, outputfile):

    if os.path.exists(outputfile):
        os.remove(outputfile)
    with open(outputfile, "w") as csvfile:
        writer = csv.writer(csvfile,dialect=csv.excel_tab)
        writer.writerow(('Total number of reads:', numb_reads))
        writer.writerow(('Total number of nucleotides:', numb_nuc))
        writer.writerow('')
        writer.writerow(('Phread quality score', 'Number of times', 'Percentage (%)'))


        od = collections.OrderedDict(sorted(dict.items()))

        for k, v in iter(od.items()):
            percentage = round((v / float(numb_nuc) * 100), 3)
            writer.writerow((k,v,percentage))

        csvfile.close()


parser = argparse.ArgumentParser(description='Script to simply analyze the qualities of a 454 file.')
parser.add_argument(dest='qualities_file', metavar='454_qualities_file', nargs=1,
            help='File with the quality values in Contigs/Scaffolds fasta file to be processed.')
parser.add_argument(dest='output_file', metavar='output_file', nargs=1, help='Output file')
args = parser.parse_args()


if __name__ == "__main__":
    print("Checking qualities of the reads.. This may take a while..")
    dict, numb_nuc, numb_reads = process_qual_file(args.qualities_file)
    print("Writing the output to " + args.output_file[0] + "..")
    write2file(dict, numb_nuc, numb_reads, args.output_file[0])
    print("Done.")
