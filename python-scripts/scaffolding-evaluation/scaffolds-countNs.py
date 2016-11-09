__author__ = 'pedro'

"""
Description: This script counts the length and number of Ns of each scaffold in a file, exporting the results in a
sorted tab separated file.
"""

#TODO import script that extracts contigs/scaffolds longer than X

import argparse
import sys
import os
from Bio import SeqIO
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
from collections import defaultdict
from collections import OrderedDict

class countN():

    def __init__(self, minlength, outputprefix, differentFasta):
        if minlength is not None:
            self.minlen = minlength
        else:
            self.minlen = 0
        self.outputprefix = outputprefix
        self.differentFasta = differentFasta
        self.finalDict = defaultdict()
        self.dictSumNs = OrderedDict()


    def processFastaFiles(self,inputFile):
        totalGenomeLength = 0
        totalGenomeLenghNoN = 0
        Ns = []
        handle = open(inputFile,'rU')
        sequences = SeqIO.parse(handle,'fasta')
        for record in sequences:

            if len(record) >= self.minlen:
                upper=record.seq.upper()
                countN = upper.count("N")
                Ns.append(countN)
                totalGenomeLength += len(record.seq)
                totalGenomeLenghNoN += len(record.seq) - countN
                if record.id in self.finalDict.keys():
                    self.finalDict[record.id].append((len(record.seq), countN,round((countN/len(record.seq))*100,2)))
                else:
                    self.finalDict[record.id] = [(len(record.seq), countN,round((countN/len(record.seq))*100,2))]
        handle.close()

        Ns.extend((totalGenomeLength,totalGenomeLenghNoN))
        file_ID = os.path.basename(inputFile).rsplit('.')[0]
        if file_ID not in self.dictSumNs:
            self.dictSumNs[file_ID] = Ns
        else:
            logging.error("%s file ID is repeated. Please check if you are not repeating the same input files.")
            exit(1)

    def processIndependentFastaFiles(self,inputFile):
        Ns = []
        total_length=0
        total_length_noN = 0
        individualFinalDict = defaultdict()
        handle = open(inputFile,'rU')
        sequences = SeqIO.parse(handle,'fasta')
        for record in sequences:

            if len(record) >= self.minlen:
                upper=record.seq.upper()
                countN = upper.count("N")
                Ns.append(countN)
                total_length += len(record.seq)
                total_length_noN += len(record.seq) - countN
                if record.id in individualFinalDict.keys():
                    logging.error("%s found more than once in %s file.." % (record.id,inputFile))
                    #individualFinalDict[record.id].append((len(record.seq), countN,round((countN/len(record.seq))*100,2)))
                else:
                    individualFinalDict[record.id] = [len(record.seq), countN,round((countN/len(record.seq))*100,2)]
        handle.close()

        outFile = os.getcwd() +'/' + str(os.path.basename(inputFile).rsplit('.')[0]) + "-indivNsInfo.txt"
        with open(outFile,'w') as outfile:
            outfile.write("#Number of Ns in file:\t%i\n" % sum(Ns))
            outfile.write("#Length of genome with Ns:\t%i\n" % total_length)
            outfile.write("#Length of genome without Ns:\t%i\n" % total_length_noN)
            outfile.write("#Fraction of Ns in genome:\t%f\n\n\n" % round((sum(Ns)/total_length)*100,2))
            outfile.write("#scaffold_id\tlength\tNs\tfraction_Ns\n")
            for k,v in sorted(individualFinalDict.values(),key = lambda x: x[0]):
                print(v)
                outfile.write(k + "\t" + '\t'.join(v))


    def writeOutput(self):
        outputfile = os.getcwd() + "/" + self.outputprefix + ".tsv"
        if os.path.exists(outputfile):
            os.remove(outputfile)
        with open(outputfile, "w") as outfile:
            outfile.write("#File\tNs in file\tlength of genome with Ns\tlength of genome without Ns\tfraction of Ns in genome\n")
            for k,v in self.dictSumNs.items():
                outfile.write(k + "\t" + sum(v[-2:]) + "\t" + v[-2] + "\t" + v[-1] + "\t" + round((sum(int(v[-2:]))/int(v[-2]))*100,2))

            outfile.write("\n\n\n#scaffold_id\t")
            for file in self.dictSumNs.keys():
                outfile.write("length\tnumber of Ns\tfraction of Ns\t")
            outfile.write("\n")
            for key,value in sorted(self.finalDict.items(), key = lambda x: x[0]):
                outfile.write(key + "\t")
                for stats in value:
                    outfile.write("\t".join(stats) + "\t")
                outfile.write("\n")
        outfile.close()

def main():

    parser = argparse.ArgumentParser(description='Script analyse the number of Ns in scaffolds. When multiple FASTA are provided, it expects the same scaffold names '
                                                 'between files (e.g check Ns distribution after some rounds of GapClosing). If each fasta represents different scaffolding'
                                                 'strategy, please set "-d" flag.')
    parser.add_argument(dest='inputFile', metavar='input_files', nargs= '+', help='Fasta/s file/s to be processed.')
    parser.add_argument('-n', '--minlength', type=int, help='Minimum length of contigs/scaffolds to process.')
    parser.add_argument('-p', '--prefix', required = True, help='[Required].Prefix name for the .tsv output file. When "-d", prefix does not take effect.')
    parser.add_argument('-d','--differentFasta', action='store_true', help = 'FASTA files are not related to each other. Individual output files per file will be generated.')
    args = parser.parse_args()

    countNobject = countN(args.minlength, args.prefix, args.differentFasta)
    for file in args.inputFile:
        logging.info("Processing %s file." % file)
        if countNobject.differentFasta:
            countNobject.processIndependentFastaFiles(file)
        else:
            countNobject.processFastaFiles(file)

    if args.differentFasta is None:
        countNobject.writeOutput()
    logging.info("Done.")

if __name__ == "__main__":
    main()


