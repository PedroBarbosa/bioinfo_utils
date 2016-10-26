import argparse
import matplotlib
matplotlib.use('Agg')
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
    Ns = []
    handle = open(inputFile,'rU')
    sequences = SeqIO.parse(handle,'fasta')
    for record in sequences:
        upper=record.seq.upper()
        countN = upper.count("N")
        Ns.append(countN)
        list_tuples.append((record.id, len(record.seq), countN))
        lenghts.append(len(record.seq))

    handle.close()
    return  list_tuples, lenghts, Ns

def writeIndividualReadLenght(tuples, inputfile):
    basename = os.path.splitext(inputfile)[0]
    outFile = basename + "-indivReadLength.txt"
    dictNsPerReadLenght = {}
    with open(outFile,'w') as out:
        out.write("#readID\tlength\tNs\n")
        for read, length, nCounts in sorted(tuples,key=itemgetter(1),reverse=True):
            out.write("%s\t%i\t%i\n" % (read,length,nCounts))
            dictNsPerReadLenght[length] = nCounts
    out.close()
    return dictNsPerReadLenght

def countsAndStats(lengths, listNs ,inputfile):
    basename = os.path.splitext(inputfile)[0]
    outFile = basename + "-statsAndCounts.txt"
    dict_counter = Counter(lengths)

    with open(outFile,'w') as out:
        ####Read lengths####
        out.write("Total number of reads:\t%i\n" % len(lengths))
        out.write("Mean read length:\t%i\n" % np.mean(lengths))
        out.write("Median of read lengths:\t%i\n" % np.median(lengths))
        out.write("Maximum read length:\t%i\n" % np.max(lengths))
        out.write("Minimum read length:\t%i\n" % np.min(lengths))
        out.write("10 most commons read lengths:\n#length\t#ocurrences\n")
        most_common = dict_counter.most_common(10)
        for val in most_common:
            out.write("%i\t%i\n" % (val[0],val[1]))


        ####Ns#####
        out.write("\n\nTotal number of Ns in reads:\t%i\n" % sum(listNs))
        out.write("Number of reads with Ns:\t%i\n" % (len(listNs) - listNs.count(0)))
        out.write("Mean number of Ns per read:\t%i\n" % np.mean(listNs))
        out.write("Median of Ns per read:\t%i\n" % np.median(listNs))
        out.write("Maximum number of Ns in one read:\t%i\n" % np.max(listNs))

        ####Counts for histogram of read lengths#####
        out.write("\nValues for histogram counts:\n")
        for key in sorted(dict_counter.keys()):
            out.write("%i\t%i\n" % (key,dict_counter[key]))

    return dict_counter



def drawHistograms(counter, dictNsReadLength, inputfile):
    basename = os.path.splitext(inputfile)[0]
    outFileReadLength = basename + "-histReadLength.png"
    outFileNsReadLength = basename + "-histNsPerReadLength.png"

    ####Read length#####
    val, weight = zip(*[(k, v) for k,v in counter.items()])
    plt.hist(val, bins= 20,weights=weight, range=[0,2000])
    plt.xlabel('Read lenght bins')
    plt.ylabel('Counts')
    plt.savefig(outFileReadLength)
    plt.close()

    ####Ns###########
    val2, weight2 = zip(*[(k, v) for k,v in dictNsReadLength.items()])
    plt.hist(val2, bins= 100,weights=weight2, range=[0,2000], color="green")
    plt.xlabel('Read lenght bins')
    plt.ylabel('Ns')
    plt.savefig(outFileNsReadLength)


parser = argparse.ArgumentParser(description='Script analyse the lenght of 454 reads. They must be converted to fasta. ')
parser.add_argument(dest='inputFile', metavar='input_files', nargs= '+', help='Fasta/s file/s to be processed.')
args = parser.parse_args()

for file in args.inputFile:
    logging.info("Processing %s file." % file)
    list_tuples, alllenghts, listNs = processFastaFiles(file)
    dict_NsReadLen = writeIndividualReadLenght(list_tuples,file)
    counter4hist =countsAndStats(alllenghts,listNs,file)
    drawHistograms(counter4hist,dict_NsReadLen,file)
    logging.info("Done.")