import argparse
import logging
import sys
import collections
import operator
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')

def createGenomeTableDict(refGenomeTable):
    gnmTbl={}
    with open(refGenomeTable, 'r') as infile:
        for line in infile:
            if len(line.split("\t")) != 2:
                print("Error, genome table files must have exactly 2 columns.")
                exit(1)
            elif line.split("\t")[0] in gnmTbl:
                print("Error, repeated scaffold id in genome table file %s" % refGenomeTable)
            else:
                gnmTbl[line.split("\t")[0]] = line.split("\t")[1].rstrip()
        #print("%i reference sequences processed" % len(gnmTbl))
    infile.close()
    return dict


def readAugustusGff(infile,outfile,dict):
    with open(infile,'r') as infl:
        with open(outfile,'w') as outfl:
            outfl.write("##gff-version 3\n##Output generated with Pedro script to convert IDs nomenclature in the 9th attributes field.\n")
            geneID=""
            transcriptID=""
            for line in infl:
                if not line.startswith("#"):
                    fields=line.split("\t")
                    if fields[2] == "gene":
                        geneID=fields[8].rstrip()
                        if fields[0] in dict:
                            outfl.write("##sequence-region " + fields[0] + " 1 " + dict[fields[0]])
                        else:
                            print("%s seq id not in genome table. Sequence region will not be present for this reference sequence." % fields[0])
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
    parser.add_argument("-g", metavar='-genomeTable', require=True,help="Genome table file to ouptut ##sequence-region in the gff. Advised by genome validation tools.")
    args = parser.parse_args()

    outdict = createGenomeTableDict(args.g)
    readAugustusGff(args.inputgff,args.outputfgff,outdict)


if __name__ == "__main__":
    main()
