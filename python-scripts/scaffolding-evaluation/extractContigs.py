import argparse
from itertools import islice
import sys
from collections import defaultdict
from os.path import expanduser

#This script extracts from a contigs/scaffolds file the ones longer than a given number

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def processFastaFiles(inputFile, n, outputFile):

        contig_seq = ""
        contig_id = ""
        keptcontigs,total_contigs = 0,0
        discardedcontigs = -1 ##to account for the first line
        with open(inputFile[0]) as file:
            for line in file:
                line=line.rstrip()
                #print contig_seq
                if line.startswith('>') and len(list(contig_seq)) >= n:

                    with open(outputFile, "a") as out:
                        out.write(contig_id + contig_seq)

                    contig_id = line
                    contig_seq = ""
                    keptcontigs+=1
                    total_contigs +=1

                elif line.startswith('>') and len(list(contig_seq)) < n:
                    contig_id = line
                    contig_seq = ""
                    discardedcontigs+=1
                    total_contigs +=1

                else:

                    contig_seq += line



            ###process last contig (after last >) ###
            if len(list(contig_seq)) > n:
                with open(outputFile, "a") as out:
                        out.write(contig_id + contig_seq)
                keptcontigs+=1
            else:

                discardedcontigs+=1

        print("\nNumber of contigs in file:\t", total_contigs)
        print("Number of contigs longer than", n , ":\t", keptcontigs)
        print("Number of contigs discarded:\t", discardedcontigs)
        print("Percentage kept:", (float(keptcontigs) / total_contigs) * 100, "%")


        file.close()
        out.close()



parser = MyParser()
parser.add_argument(dest='inputFile', metavar='input_files', nargs=1,
                    help='Contig Fasta file to be processed.')
parser.add_argument('-n', '--contiglength', metavar='', type=int, required=True, help='Minimum length of contigs to extract')
parser.add_argument('-o', '--outputFile', metavar='', required = True, help='Output file')
args = parser.parse_args()


processFastaFiles(args.inputFile,args.contiglength,args.outputFile)

