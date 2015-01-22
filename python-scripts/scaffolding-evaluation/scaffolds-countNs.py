__author__ = 'pedro'

"""
Description: This script counts the length and number of Ns of each scaffold in a file, exporting the results in a
sorted tab separated file.
"""

#TODO import script that extracts contigs/scaffolds longer than X

import argparse
import sys
import os
import csv
#from extractContigs import processFastaFiles


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        self.print_help()
        sys.exit(2)

#class countN():
#    def __init__(self,scaffolds, minlength, outputfile):
#        self.scaff = scaffolds
#        self.minlen = minlength
#        self.output = outputfile

def count(inputfile):

        scaffold_id = ""
        map={}
        n,total = 0,0
        firstScaff=True
        with open(inputfile[0]) as file:
            for line in file:
                if line.startswith('>'):
                    if firstScaff:
                        scaffold_id = line
                        firstScaff = False
                    else:

                        map[scaffold_id] = (total,n)
                        #print(contig_id,n,total)
                        scaffold_id = line
                        n,total = 0,0
                else:
                     n += line.count("N")
                     total += len(line)

            ###process last scaffold (after last >)
            map[scaffold_id] = (total,n)



            file.close()
            return map


def write2file(dict,outputfile):

    if os.path.exists(outputfile):
        os.remove(outputfile)
    with open(outputfile, "w") as csvfile:
        writer = csv.writer(csvfile,dialect=csv.excel_tab)
        writer.writerow(('Scaffolds_id','Length','Number of Ns','Proportion of Ns (%)'))

        for n, total in sorted(dict.values(), key=lambda e: e[0], reverse=True):
            tuple = (n, total)

            #get keys of the tuples sorted by the length of scaffolds and write to file
            for key in dict.keys():
                if dict[key] == tuple :
                    key.strip()
                    percentageNs = round((float(tuple[1]) / tuple[0]) * 100,2)
                    writer.writerow((key.rstrip(), str(tuple[0]),str(tuple[1]), str(percentageNs)))

                else:
                    continue

    csvfile.close()


parser = MyParser()
parser.add_argument(dest='inputFile', metavar='input_files', nargs=1,
            help='Contigs/Scaffolds fasta file to be processed.')
parser.add_argument('-n', '--minlength', metavar='', type=int, required=False, help='Minimum length of contigs/scaffolds to process')
parser.add_argument('-p', '--prefix', metavar='', required = False, help='Prefix name for the .csv output file')
args = parser.parse_args()



if __name__ == "__main__":
    print("Counting ocurrences in scaffolds..")
    dict = count(args.inputFile)
    print("Writing output to csv..")
    #need to be .tsv file in order for windows and unix interpret it correctly
    if args.prefix:
        write2file(dict, args.prefix + ".tsv")
    else:
        write2file(dict,"listNs.tsv")

