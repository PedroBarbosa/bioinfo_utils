import argparse
import logging
import sys
import collections
import operator
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')


def readAugustusGff(infile,outfile):
    with open(infile,'r') as infl:
        with open(outfile,'w') as outfl:
            outfl.write("#Output generated with Pedro script to convert IDs nomenclature in the 9th attributes field.\n")
            geneID=""
            transcriptID=""
            for line in infl:
                if not line.startswith("#"):
                    fields=line.split("\t")
                    if fields[2] == "gene":
                        geneID=fields[8].rstrip()
                        outfl.write('\t'.join(fields[:-1]) + '\tID=' + geneID + "\n")
                    elif fields[2] == "transcript":
                        transcriptID=fields[8].rstrip()
                        outfl.write('\t'.join(fields[:-1]) + '\tID=' + transcriptID + ";Parent=" + geneID + "\n")
                    elif fields[2] == "start_codon" or fields[2] == "stop_codon" or fields[2] == "exon" or fields[2] == "CDS":
                         outfl.write('\t'.join(fields[:-1]) + '\tParent=' + transcriptID + "\n")
                    elif fields[2] == "intron":
                        continue

    infl.close()
    outfl.close()



def main():
    parser = argparse.ArgumentParser(description='Script to convert native non gff3 augustus format into a readable gff3.')
    parser.add_argument(dest='inputgff', metavar='augustusGff', help='Input file.')
    parser.add_argument(dest='outputfgff', metavar='outgff3', help='Name of the output file.')
    args = parser.parse_args()
    readAugustusGff(args.inputgff,args.outputfgff)


if __name__ == "__main__":
    main()
