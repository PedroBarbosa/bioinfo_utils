import argparse
import matplotlib.pyplot as plt
from collections import Counter
from operator import itemgetter
from Bio import SeqIO
import logging
import sys
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
import os
import numpy as np

def processFastaFiles(inputFile):

    list_tuples = []
    lenghts = []
    handle = open(inputFile,'rU')
    sequences = SeqIO.parse(handle,'fasta')
    for record in sequences:
        list_tuples.append((record.id, len(record.seq)))
        lenghts.append(len(record.seq))

    handle.close()
    return  list_tuples, lenghts

def writeIndividualReadLenght(tuples, inputfile):
    basename = os.path.splitext(inputfile)[0]
    outFile = basename + "-indivReadLength.txt"
    with open(outFile,'w') as out:
        out.write("#readID\tlength\n")
        for read, length in sorted(tuples,key=itemgetter(1),reverse=True):
            out.write("%s\t%i\n" % (read,length))
    out.close()


def countsAndStats(lengths, inputfile):
    basename = os.path.splitext(inputfile)[0]
    outFile = basename + "-statsAndCounts.txt"
    dict_counter = Counter(lengths)

    with open(outFile,'w') as out:
        out.write("Total number of reads:\t%i\n" % len(lengths))
        out.write("Mean read length:\t%i\n" % np.mean(lengths))
        out.write("Median of read lengths:\t%i\n" % np.median(lengths))
        out.write("Maximum read length:\t%i\n" % np.max(lengths))
        out.write("Minimum read length:\t%i\n" % np.min(lengths))
        out.write("10 most commons read lengths:\n#length\t#ocurrences\n")
        most_common = dict_counter.most_common(10)
        for val in most_common:
            out.write("%i\t%i\n" % (val[0],val[1]))

        out.write("\nValues for histogram counts:\n")
        for key in sorted(dict_counter.keys()):
            out.write("%i\t%i\n" % (key,dict_counter[key]))

    return dict_counter



def drawHistogram(counter, inputfile):
    basename = os.path.splitext(inputfile)[0]
    outFile = basename + "-histogramReadLength.png"

    val, weight = zip(*[(k, v) for k,v in counter.items()])
    plt.hist(val, bins= 15,weights=weight, range=[0,2000])
    plt.xlabel('Read lenght bins')
    plt.ylabel('Counts')
    plt.savefig(outFile)


parser = argparse.ArgumentParser(description='Script analyse the lenght of 454 reads. They must be converted to fasta. ')
parser.add_argument(dest='inputFile', metavar='input_files', nargs= '+', help='Fasta/s file/s to be processed.')
args = parser.parse_args()

for file in args.inputFile:
    logging.info("Processing %s file." % file)
    list_tuples, alllenghts = processFastaFiles(file)
    writeIndividualReadLenght(list_tuples,file)
    counter4hist =countsAndStats(alllenghts,file)
    drawHistogram(counter4hist,file)