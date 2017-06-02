from Bio import SeqIO
import argparse
import sys

def createDict(auxFile,column):
    maps={}
    with open(auxFile,'r') as infile:
        for l in infile:
            line=l.rstrip()
            if not line.startswith("#") and len(line.split("\t")) != 2:
                print("Two column mapping file is required.Exiting. ")
                exit(1)
            else:
                if line.split("\t")[column - 1] in maps:
                    print("Error. %s read header is duplicate.")
                    exit(1)
                else:
                    if column == 1:
                        maps[line.split("\t")[column-1]] = line.split("\t")[column]
                    elif column == 2:
                        maps[line.split("\t")[column-1]] = line.split("\t")[column-2]

    infile.close()
    return maps

def replaceFastqHeaders(fastq,mappingDict):
    handle = open(fastq,"rU")
    for seq_record in SeqIO.parse(handle,"fastq"):
        old_header = seq_record.id
        print(old_header)

def main():
    parser = argparse.ArgumentParser(description='Script to rename the header of a Fa file, either by adding string in the end of header or totally replace it. Biopython needed')
    parser.add_argument(dest='fastqFile', metavar='fastqFile', help='FASTQ file to process.')
    parser.add_argument('-a', metavar='auxiliarMapping', required=True, help='Auxiliar 2 columns tab file mapping read names to new IDs.')
    parser.add_argument('-c', metavar='columnOfOldID',type=int,default=2, choices=(1,2),help='Column where old fastq headers are displayed. Default: Second, represented by 2.')
    args = parser.parse_args()

    mapsDict = createDict(args.auxiliarMapping,args.columnOfNewID)
    replaceFastqHeaders(args.fastqFile,mapsDict)
