from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import argparse
import sys

def createDict(auxFile,column):
    print("Hashing mapping file..")
    maps={}
    with open(auxFile,'r') as infile:
        for l in infile:
            line=l.rstrip()
            if not line.startswith("#") and len(line.split("\t")) != 2:
                print("Two column mapping file is required.Exiting.")
                exit(1)
            else:
                if line.split("\t")[column - 1] in maps:
                    print("Error. %s read header is duplicate.")
                    exit(1)
                else:
                    if column == 1:
                        maps[line.split("\t")[column-1]] = line.split("\t")[column]
                    elif column == 2:
                        print(line.split("\t")[column-1])
                        maps[line.split("\t")[column-1]] = line.split("\t")[column-2]
    
    print("%i read names hashed" % len(maps))
    infile.close()
    return maps

def replaceFastqHeaders(fastq,mappingDict):
    print("Processing fastq file..")
    handle = open(fastq,"rU")
    handle_out=open("out.fastq","w")
    extraReads=0
    for (title, sequence, quality) in FastqGeneralIterator(handle):
        if title not in mappingDict:
            extraReads+=1
        else:
            handle_out.write("@%s\n%s\n+\n%s\n" % (mappingDict[title], sequence, quality))
    print("Number of reads absent from the mapping file, thus not reported in the output file:\t%i" % extraReads)

def replaceUpToFirstSpace(fastq):
    print("Processing fastq file..")
    handle=open(fastq,"rU")
    handle_out=open("out.fastq","w")
    for (title,sequence,quality) in FastqGeneralIterator(handle):
        newHeader = title.split(" ")[0]
        handle_out.write("@%s\n%s\n+\n%s\n" % (newHeader, sequence, quality))
    handle.close()
    handle_out.close()

def main():
    parser = argparse.ArgumentParser(description='Script to rename the header of a Fa file, either by adding string in the end of header or totally replace it. Biopython needed')
    parser.add_argument(dest='fastqFile', metavar='fastqFile', help='FASTQ file to process.')
    parser.add_argument('-w','--removeSpace',action='store_true',help='Flag to select first string up to the first whitespace.')
    parser.add_argument('-a', metavar='auxiliarMapping', help='Use of an auxiliar 2 columns tab file to replace old ids with the new ones.')
    parser.add_argument('-c', metavar='columnOfOldID',type=int,default=2, choices=(1,2),help='When using -a argument, this refers to the column where old fastq headers are displayed. Default: Second, represented by 2.')
    args = parser.parse_args()
    
    if args.a and not args.removeSpace:
        mapsDict = createDict(args.a,args.c)
        replaceFastqHeaders(args.fastqFile,mapsDict)
    elif args.removeSpace:
        replaceUpToFirstSpace(args.fastqFile)

if __name__ == "__main__":
    main()
