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
        list_tuples.append((record.id, len(record.seq),countN))
        lenghts.append(len(record.seq))

    handle.close()
    return  list_tuples, lenghts, Ns

def processQualityFiles(qualFile):

    handle=open(qualFile.rstrip(),'rU')
    qualities = SeqIO.parse(handle,'qual')
    dict_qualities = {}
    for record in qualities:
        read_qualities = record.letter_annotations["phred_quality"]
        dict_qualities[record.id] = np.mean(read_qualities)

    return  dict_qualities

def writeIndividualReadLenght(tuples, dict_qual, inputfile):
    basename = os.path.basename(inputfile).rsplit('.',1)[0]
    outFile = os.getcwd() +'/' + str(basename) + "-indivReadInfo.txt"

    dictNsPerReadLenght = {}
    dictFor4QualScatterPlot = {}
    with open(outFile,'w') as out:
        out.write("#readID\tlength\tNs\t%Ns\tAvgPhreadQuality\n")
        for read, length, nCounts in sorted(tuples,key=itemgetter(1),reverse=True):
            out.write("%s\t%i\t%i\t%f\t%f\n" % (read,length,nCounts,round((nCounts/length)*100,4),dict_qual[read]))
            dictNsPerReadLenght[length] = nCounts
            dictFor4QualScatterPlot[length] = dict_qual[read]
    out.close()
    return dictNsPerReadLenght, dictFor4QualScatterPlot

def countsAndStats(lengths, listNs ,dict_qual, inputfile):
    basename = os.path.basename(inputfile).rsplit('.',1)[0]
    outFile = os.getcwd() +'/' + str(basename) + "-statsAndCounts.txt"
    dict_counter = Counter(lengths)

    with open(outFile,'w') as out:
        ####Read lengths####
        out.write("Total number of reads:\t%i\n" % len(lengths))
        out.write("Mean read length:\t%f\n" % np.mean(lengths))
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
        out.write("Mean number of Ns per read:\t%f\n" % np.mean(listNs))
        out.write("Median of Ns per read:\t%i\n" % np.median(listNs))
        out.write("Maximum number of Ns in one read:\t%i\n" % np.max(listNs))


        ####Qualities#####
        out.write("\nAverage read phread quality:\t%f" % np.mean(list(dict_qual.values())))
        out.write("\nMedian of read phread quality:\t%i" % np.median(list(dict_qual.values())))
        out.write("\nLowest average quality of a read:\t%s" % np.min(list(dict_qual.values())))
        out.write("\nHighest average quality of a read:\t%s" % np.max(list(dict_qual.values())))

        ####Counts for histogram of read lengths#####
        out.write("\n\n\nValues for histogram counts:\n")
        for key in sorted(dict_counter.keys()):
            out.write("%i\t%i\n" % (key,dict_counter[key]))

    return dict_counter



def drawHistograms(counter, dictNsReadLength, inputfile):
    basename = os.path.basename(inputfile).rsplit('.',1)[0]
    outFileReadLength = os.getcwd() +'/' + str(basename) + "-histReadLength.png"
    outFileNsReadLength = os.getcwd() +'/' + str(basename) + "-histNsPerReadLength.png"

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
    plt.close()

def drawScatterPlot(dictFor4QualScatterPlot,inputfile):
    basename = os.path.basename(inputfile).rsplit('.',1)[0]
    outFileScatterQual = os.getcwd() +'/' + str(basename) + "-scatterQualities.png"

    for length,qual in dictFor4QualScatterPlot.items():
        plt.scatter(length,qual,color='grey',s=25,alpha=0.6)
    plt.xlabel("Read lenght")
    plt.ylabel("Average Phread Qualities")
    plt.xlim(0,2000)
    plt.savefig(outFileScatterQual)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Script analyse the lenght of 454 reads. They must be converted to fasta. ')
    parser.add_argument(dest='inputFile', metavar='input_files', nargs= '+', help='Fasta/s file/s to be processed.')
    parser.add_argument('-l', '--listqualityFiles', required=True, help='Add quality analysis for each fasta. Files must have the same name as the fasta, except for the extension.')
    args = parser.parse_args()

    for file in args.inputFile:
        logging.info("Processing %s file." % file)
        if args.listqualityFiles:

            qualities=""
            with open(args.listqualityFiles,'r') as listFiles:
                for qualFile in listFiles:
                    if os.path.abspath(file).rsplit('.', 1)[0] == os.path.abspath(qualFile).rsplit('.',1)[0]:
                        qualities = qualFile.rstrip()
                        break
            if qualities=="":
                logging.error("Error. %s file does not have quality file associated. Please check carefully the name of the files." % file)
                logging.info("Problematic FASTA file path without extension:\t%s" % os.path.abspath(file).rsplit('.', 1)[0])
                exit(1)

            logging.info("Taking a look at quality file (%s)." % qualities)
            dict_qual = processQualityFiles(qualities)

        logging.info("Checking quality..Done.")
        list_tuples, alllenghts, listNs = processFastaFiles(file)
        logging.info("Calculating overall stats..")
        dict_NsReadLen,dictFor4QualScatterPlot = writeIndividualReadLenght(list_tuples,dict_qual,file)
        counter4hist =countsAndStats(alllenghts,listNs,dict_qual,file)
        drawHistograms(counter4hist,dict_NsReadLen,file)
        drawScatterPlot(dictFor4QualScatterPlot, file)
        logging.info("Done.")


if __name__ == "__main__":
    main()