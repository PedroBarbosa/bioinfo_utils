#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import SeqIO
import argparse

def processLastal(lastOut):
    with open(lastOut,"r") as infile:

        problematicList=set()
        pairwiseAln =set()
        for l in infile:
            line=l.rstrip()
            if not line.startswith("#") and len(line.split()) > 10:

                (score, query, qstart, qalgsize, qstrand, qsize, reference, refstart, refalgsize, rstrand, refsize, blocks) = line.split()[:12]
                (score, qstart, qalgsize, qsize, refstart, refalgsize, refsize) = map(int, (score, qstart, qalgsize, qsize, refstart, refalgsize, refsize))


                #removed score as it seems even in perfect matches the score is not always the same
                if all(v == qalgsize for v in (qalgsize, qsize, refalgsize, refsize)) and query != reference:
                    if query.isdigit():
                        pairwiseAln.add((query,reference))
                    elif reference.isdigit():
                        pairwiseAln.add((reference,query))
                    else:
                        print("None of fasta header seems to represent redundans output (integer values)")
                        exit(1)

                else:
                    problematicList.add(query)
                    #if refsize/qsize > 0.8:
                    #    print(qsize,refsize)

                    #    print(qsize)
                    #    print(qalgsize,refalgsize)
                    #    print(str(score) + "\n")
                    #print(line)
                    #print(qalgsize,qsize,refalgsize,refsize)
    infile.close()
    outdict = dict((x, y) for x, y in pairwiseAln)
    print("Number of perfect matches:\t%i" % len(outdict))
    
    listq=set()
    for id in problematicList:
        if not id in outdict:
            listq.add(id)
    return outdict, listq

def writeFileProblematic(listq):
    outfile=open("problematicQuery.txt", 'w')
    outfile.write('\n'.join(listq))

def processFasta(dicPairs, redundands,outputfile):

    handle=open(redundands, "rU")
    redundands_seq = SeqIO.parse(handle,'fasta')
    output = open(outputfile,'w')
    for record in redundands_seq:
        if record.id in dicPairs:
            output.write(">" + dicPairs[record.id] + "\n" + str(record.seq) + "\n")
        else:
            print("Redundans id %s not present in pairwise dictionary." % record.id)
            #exit(1)

    handle.close()
    output.close()


def main():
    parser = argparse.ArgumentParser(description='Script to replace fasta headers generated by redundans by the previous ones.')
    parser.add_argument(dest='lastOutput', metavar='lastTabOut', help='Lastal output file in TAB format mapping the reference inputted to redundands against redundans output.')
    parser.add_argument(dest='fasta_file', metavar='fastaFile', help='Redundans output file.')
    parser.add_argument(dest='outputFile', metavar='outputFile', help='File to write corrected output.')
    args = parser.parse_args()

    print("Processing lastal output..")
    validPairWise,listq = processLastal(args.lastOutput)
    writeFileProblematic(listq)
    print("Writing new fasta..")
    processFasta(validPairWise,args.fasta_file,args.outputFile)

if __name__ == "__main__":
    main()
