import argparse
from itertools import islice
import sys
import re
import math
from os.path import expanduser

#This script predicts which ORFs may represent full length genes and writes two different files based on the complete and incomplete putative ORFs

__author__ = 'pedro'

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def calculateFullGenes(inputFile, outputDirectory):
    print "Processing nucleotide ORFs file ..."
    isFullGene = False
    firstLine = True
    gene_id, seq = "", ""
    full, incomplete = 0 ,0
    with open(inputFile[0]) as file:
            for line in file:

                if line.startswith('>') and firstLine == False:

                    seq.strip()
                    initiation = seq[:3]
                    stop = seq[-4:]

                    stop_real = re.sub("[^a-z0-9]+","", stop, flags=re.IGNORECASE)

                    if initiation.upper() == "ATG" and stop_real.upper() in set(['TAG','TAA','TGA']):
                        isFullGene = True

                    else:
                        isFullGene = False


                    if isFullGene:
                        with open(outputDirectory + "/complete_ORFs.fasta", "a") as out1:
                            out1.write(gene_id + seq)
                        out1.close()
                        full+=1

                    else:
                        with open(outputDirectory + "/incomplete_ORFs.fasta", "a") as out2:
                            out2.write(gene_id + seq)
                        out2.close()
                        incomplete+=1


                    gene_id = line
                    seq = ""

                elif firstLine == False:
                    seq += line


                else:
                    gene_id = line
                    firstLine = False


            ###process last ORF (after last >) ###
            initiation = seq[:3]
            stop = seq[-4:]
            stop_real = re.sub("[^a-z0-9]+","", stop, flags=re.IGNORECASE)
            if initiation.upper() == "ATG" and stop_real.upper() in set(['TAG','TAA','TGA']):
                isFullGene = True
            else:
                isFullGene = False

            if isFullGene:
                with open(outputDirectory + "/complete_ORFs.fasta", "a") as out1:
                    out1.write(gene_id + seq)
                out1.close()
                full+=1

            else:
                with open(outputDirectory + "/incomplete_ORFs.fasta", "a") as out2:
                    out2.write(gene_id + seq)
                out2.close()
                incomplete+=1


    print "\nTotal number of ORFs in the file:\t" , full + incomplete
    print "Number of full length ORFs predicted:\t" ,  full
    print "Number of incomplete ORFs predicted:\t" , incomplete
    print "Percentage of full length ORFs:\t", (float(full)/(full+incomplete)) * 100, "%"



parser = MyParser()
parser.add_argument(dest='inputFile', metavar='input_file', nargs=1, help='File with the predicted ORFs in the nucleotide alphabet')
parser.add_argument('-o', '--outputDirectory', metavar='', default = expanduser("~"), help='Output directory of the files. Default: $HOME')
args = parser.parse_args()

if __name__ == "__main__":



    ####check file alphabet #############
    with open(args.inputFile[0], "r") as file:
        next_2_lines = list(islice(file, 2))
	m = re.search('(^[agtcAGTC]+$)',next_2_lines[1])
        if not m:
            sys.stderr.write('\nIs your file with the DNA alphabet ? check it out!\nExiting..\n')
            sys.exit(2)

    calculateFullGenes(args.inputFile, args.outputDirectory)
