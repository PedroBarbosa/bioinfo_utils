import argparse
from lib2to3.pgen2.tokenize import double3prog
from decimal import *

__author__ = 'pedro'

parser = argparse.ArgumentParser(description='This is a script to correct the evalue of rapsearch results')
parser.add_argument('-i', dest='blast_results_file', help='Blast results file', required=True)
parser.add_argument('-e', dest='evalue_threshold', help='Threshold to accept hits', required=True)
parser.add_argument('-o', dest='outputed_file', help='Corrected evalue file', required=True)
args = parser.parse_args()

def processEvalue(blastFile, eval_threshold):
    print("Reading and processing file ...")
    filename = args.outputed_file
    open(filename, 'w').close() ## empty output file, just in case
    fileoutput = open (filename, 'w')
    with open(blastFile) as file:
        for line in file:
            if not line.startswith('#'):
                vector = line.split("\t") #read itself
                eValue = round(10 ** (Decimal(vector[10])),8) #new eValue
                if eValue <= float(eval_threshold):
                    vector[10] = str(eValue)
                    fileoutput.write('\t'.join(vector))
            else:
                vector = line.split("\t") #read itself
                if len(vector) > 10:
                    vector[10] = "e-value"
                    fileoutput.write('\t'.join(vector))
                else:
                    fileoutput.write(line)
    print("Done.")
    return fileoutput

processEvalue(args.blast_results_file, args.evalue_threshold)
