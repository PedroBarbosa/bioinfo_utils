import argparse
import matplotlib
matplotlib.use('Agg')
from Bio import SeqIO
from collections import OrderedDict
import os

def declareGLobal():
    global dataStruct
    dataStruct=OrderedDict(dict)

def processFastaFiles(inputFile):

    smrtcell_id = os.path.basename(inputFile).split("_")[3]

    if smrtcell_id in dataStruct:
        print("Processing another file of existing %s smrtcell" % smrtcell_id)
        indivSmrtcellDic = dataStruct[smrtcell_id]
    else:
        print("Processing first file of %s smrtcell" % smrtcell_id)
        indivSmrtcellDic = OrderedDict(list)
    rq = 0

    handle = open(inputFile,'rU')
    sequences = SeqIO.parse(handle,'fasta')
    for record in sequences:

        hole_number=record.id.split("/")[2]
        if not "RQ=" in record.id:
            print("read quality tag may not be present in fasta headers.")
        else:
            rq=record.id.split("RQ=")[1]

        if not hole_number in indivSmrtcellDic:
            indivSmrtcellDic[hole_number] = [rq]

        indivSmrtcellDic[hole_number].append(len(record.seq))


    dataStruct[smrtcell_id] = indivSmrtcellDic
    handle.close()


def main():

    parser = argparse.ArgumentParser(description='Script to check overall stats from pacbio subreads data in fasta format.')
    parser.add_argument(dest='fasta_file', metavar='fasta', nargs="+", help='Fasta files to be processed.')
    parser.add_argument('-o', metavar = 'output_basename',required=True,help='Basename to write the output files.')
    parser.add_argument('-l', '--list', action='store_true', help='Input is a file listing all the fasta together, one per line.')
    args = parser.parse_args()

    declareGLobal()
    if len(args.fasta_file) > 1 and args.list:
        print("When '-l' set, please don't provide more than 1 file in the positional arguments." )
        exit(1)

    elif args.list:
        with open(args.fasta_file[0], 'r') as listFasta:
            for line in listFasta:
                processFastaFiles(line)

    else:
        for file in args.fasta_file:
            processFastaFiles(file)

if __name__ == "__main__":
    main()