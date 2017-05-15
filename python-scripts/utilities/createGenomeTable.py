import argparse
import os
import collections
from collections import OrderedDict
from Bio import SeqIO
def processFastaFiles(inputFile,removeSpaces):

    contig_seq = ""
    contig_id = ""
    final_dict = {}
    handle=open(inputFile, "rU")
    sequences = SeqIO.parse(handle,'fasta')
    for record in sequences:
        if removeSpaces:
            final_dict[record.description.split()[0]] = len(record.seq)
        else:
            final_dict[record.description] = len(record.seq)

    handle.close()
    return  final_dict


def writeOutput(final_dict,orderedFasta,outputFile, inputFileName):

    sorted_dict = OrderedDict(sorted(final_dict.items(),key = lambda x: x[1][0],reverse=True))
    #print(type(sorted_dict))
    with open(outputFile, "w") as out:
        #for k,v in sorted(final_dict.items(),key = lambda x: x[1][0],reverse=True):
        for k,v in iter(sorted_dict.items()):
            out.write(k + "\t" + str(v[0]) + "\n")

    if orderedFasta:
        outFastaFile = os.getcwd() +'/' + str(os.path.basename(inputFileName).rsplit('.')[0]) + "-ordered.fasta"
        with open(outFastaFile, "w") as outFasta:
            for k,v in iter(sorted_dict.items()):
                outFasta.write('>' + k + "\n" + str(v[1]) + "\n")


parser = argparse.ArgumentParser(description='Script to create genome table from fasta file. [chromossome_id length_bp].')
parser.add_argument(dest='inputFile', metavar='input_files', nargs=1,
                    help='Fasta file to be processed.')
parser.add_argument('-o', '--outputFile', metavar='', required = True, help='Output file')
parser.add_argument('-f', '--fastaOrdered',action='store_true', help='Flag that writes new fasta ordered by legnth of sequences.' )
parser.add_argument('-r', '--removeSpaces',action='store_true', help='Flag that writes output file with headers truncated up to the first space found.' )
args = parser.parse_args()


dict_out = processFastaFiles(args.inputFile[0], args.removeSpaces)
writeOutput(dict_out,args.fastaOrdered,args.outputFile,args.inputFile[0])