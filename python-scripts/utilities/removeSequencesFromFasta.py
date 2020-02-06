#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import SeqIO
import argparse

def processFasta(fastafile,removeList,outputfile):
    toRemove = set()
    with open(removeList) as f:
        for line in f:
            line = line.strip()
            if line != "" and ">" in line:
                toRemove.add(line[1:])
            elif line != "":
                toRemove.add(line)

    handle=open(fastafile, "rU")
    fasta_sequences = SeqIO.parse(handle,'fasta')
    output_handle = open(outputfile, "w")
    for record in fasta_sequences:
        if record.id not in toRemove and len(record.seq) > 0:
            SeqIO.write(record,output_handle, "fasta")


def main():
    parser = argparse.ArgumentParser(description='Script to remove some sequences from FASTA file based on a list of headers, one per line.')
    parser.add_argument(dest='fasta_file', metavar='fastaFile', help='Fasta file file to process.')
    parser.add_argument(dest='listHeaders', metavar='headersFile', help='File with fasta headers to remove.')
    parser.add_argument(dest='outputFile', metavar='outputFile', help='File to write output.')
    args = parser.parse_args()

    print(args)
    processFasta(args.fasta_file,args.listHeaders,args.outputFile)

if __name__ == "__main__":
    main()