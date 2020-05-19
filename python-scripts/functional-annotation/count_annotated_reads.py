import argparse
from collections import defaultdict
import subprocess
import os

__author__ = 'pedro'


parser = argparse.ArgumentParser(description='This is a script to count the reads annotated in a sample')
parser.add_argument('-i', dest='blast_results_file', help='Blast results file', required=True)
parser.add_argument('-r', dest='reads_file', help='Original reads file', required=True)
parser.add_argument('-f', dest='reads_format', required=True, default= 'fasta', help = 'Format of the input reads file: fasta or fastq. Default: fasta')
parser.add_argument('-g', dest='blastFile_format', required=True, default= 'blastTAB', help = 'Format of the BLAST file: blastTAB, blastX, blastP, rpsblast or maf. Default: blastTAB')

args = parser.parse_args()

def countReadsFromFasta(fastaReadsFile):
    print("Counting number of reads ...")
    ps = subprocess.Popen(['grep', '>', fastaReadsFile], stdout=subprocess.PIPE)
    return subprocess.check_output(['wc', '-l'], stdin=ps.stdout,universal_newlines=True)
    ps.wait()
    #return subprocess.Popen(["head", readsFile], stdout=subprocess.PIPE).communicate()[0]

def countReadsFromFastq(fastqReadsFile):
    print("Counting number of reads ...")
    ps = subprocess.Popen(['cat', fastqReadsFile],stdout=subprocess.PIPE)
    reads = subprocess.check_output(['wc', '-l'], stdin=ps.stdout,universal_newlines=True)
    ps.wait()
    return int(reads) / 4

def countAnnotatedReadsBlastTab(blastFile):
    print("Counting number of annotated reads ...")
    annotated_reads = 0
    previous_line = ""
    with open(blastFile) as file:
        for line in file:

            if not line.startswith('#'):
                read = line.split("\t")[0] #read itself
                if read == previous_line: #if is different read
                    previous_line = read
                else:
                    annotated_reads += 1
                    previous_line = read

    return annotated_reads

def countAnnotatedReadsBlastX(blastFile):
    print("Counting number of annotated reads ...")
    annotated_reads = 0
    previous_line = ""
    with open(blastFile) as file:
        for line in file:

            if line.startswith('Query='):
                  annotated_reads += 1

    return annotated_reads

def countAnnotatedReadsBlastP(blastFile):
    print("Counting number of annotated reads ...")
    annotated_reads = 0
    hasHit = True
    firstgene = True
    with open(blastFile) as file:
        for line in file:

            if line.startswith('Query=') and firstgene == True:

                firstgene = False

            elif line.startswith('Query=') and firstgene == False:

                if hasHit:
                    annotated_reads += 1
                else:
                    hasHit = True

            elif "No hits" in line:

                hasHit = False

            else:
                continue

        ##process last gene
        if hasHit:
            annotated_reads +=1

    return annotated_reads

def countAnnotatedReadsRPSBlast(blastfile):
    print("Counting number of annotated reads ...")
    annotated_reads = 0
    firstLine = True
    hasHit = False

    with open(blastfile) as file:
        for line in file:
             if line.startswith('Query=') and firstLine == False:

                    if hasHit:
                        annotated_reads += 1

                    else:
                        hasHit = True

             elif "No hits" in line:

                    hasHit = False

             elif firstLine == True:
                 firstLine = False

             else:
                 continue

        ##process last gene
        if hasHit:
            annotated_reads +=1


    return annotated_reads


def countAnnotatedReadsMaf(blastFile):
    print("Counting number of annotated reads ...")
    annotated_reads = 0
    isFirst_s = True
    dict = defaultdict(int)
    with open(blastFile) as file:
        for line in file:
            if line.startswith('s') and isFirst_s:
                isFirst_s = False
            elif line.startswith('s'):
                #print line
                read = line.split(" ")[1]
                #print read
                isFirst_s = True
                if read not in dict:
                    dict[read] = 1
                else:
                    dict[read] = dict.get(read) + 1

    annotated_reads = len(dict)
    return annotated_reads


def calculatePercentage(totalReads,annotatedReads):
    percentage = round(float(annotatedReads) / float(totalReads) * 100,4)
    print("Number of reads: ", totalReads.strip())
    print("Reads with blast hits: ", annotatedReads)
    print("Percentage of the reads annotated in the sample: ", percentage, "%")
    filename = args.blast_results_file + "_stats.txt"
    fileoutput = open (filename, 'w')
    fileoutput.write("Number of reads:\t" + totalReads.strip() + "\n")
    fileoutput.write("Reads with blast hits:\t" +  str(annotatedReads) +"\n")
    fileoutput.write("Percentage of the reads annotated in the sample:\t" + str(percentage) +"%\n")
    fileoutput.close()


if args.reads_format == "fasta":
    if args.blastFile_format == "blastTAB":
        reads = countReadsFromFasta(args.reads_file)
        annotated = countAnnotatedReadsBlastTab(args.blast_results_file)
        calculatePercentage(reads,annotated)
    elif args.blastFile_format == "blastX":
        reads = countReadsFromFasta(args.reads_file)
        annotated = countAnnotatedReadsBlastX(args.blast_results_file)
        calculatePercentage(reads,annotated)
    elif args.blastFile_format == "rpsblast":
        reads = countReadsFromFasta(args.reads_file)
        annotated = countAnnotatedReadsRPSBlast(args.blast_results_file)
        calculatePercentage(reads,annotated)
    elif args.blastFile_format == "blastP":
        reads = countReadsFromFasta(args.reads_file)
        annotated = countAnnotatedReadsBlastP(args.blast_results_file)
        calculatePercentage(reads,annotated)
    elif args.blastFile_format == "maf":
        reads = countReadsFromFasta(args.reads_file)
        annotated = countAnnotatedReadsMaf(args.blast_results_file)
        calculatePercentage(reads,annotated)
    else:
        raise argparse.ArgumentTypeError("Format of the BLAST file not known")



elif args.reads_format == "fastq":
    if args.blastFile_format == "blastTAB":
        reads = countReadsFromFastq(args.reads_file)
        annotated = countAnnotatedReadsBlastTab(args.blast_results_file)
        calculatePercentage(reads,annotated)
    elif args.blastFile_format == "blastX":
        reads = countReadsFromFastq(args.reads_file)
        annotated = countAnnotatedReadsBlastX(args.blast_results_file)
        calculatePercentage(reads,annotated)
    elif args.blastFile_format == "rpsblast":
        reads = countReadsFromFastq(args.reads_file)
        annotated = countAnnotatedReadsRPSBlast(args.blast_results_file)
        calculatePercentage(reads,annotated)
    elif args.blastFile_format == "blastP":
        reads = countReadsFromFastq(args.reads_file)
        annotated = countAnnotatedReadsBlastP(args.blast_results_file)
    elif args.blastFile_format == "maf":
        reads = countReadsFromFastq(args.reads_file)
        annotated = countAnnotatedReadsMaf(args.blast_results_file)
        calculatePercentage(reads,annotated)
    else:
        raise argparse.ArgumentTypeError("Format of the BLAST file not known")

else:
    raise argparse.ArgumentTypeError('Format of the reads file not known')
