import argparse
from Bio import SeqIO

def processAuxFile(auxfile):
    #print("Reading auxiliar mapping file..")
    dict={}
    with open(auxfile,'r') as infile:
        for line in infile:
            if len(line.split("\t")) == 2:
                locus=line.split("\t")
                if not locus[0] in dict.keys():
                    dict[line.split("\t")[0]] = line.split("\t")[1].rstrip()
                else:
                    print("Repeated locus names in auxiliar files. Please take a look on this.")
                    exit(1)
            else:
                print("Problems with input auxiliar file. Please provide 2 columns in this file.")
                exit(1)
    return dict

def processGenbank(gbk,out,dict):

    #print("Replacing genbank records..")
    with open(out,'w') as outfile:
        for gb_record in SeqIO.parse(open(gbk,"r"), "genbank"):
            gb_record.name = dict[gb_record.name]
            SeqIO.write(gb_record, outfile, "genbank")
    #print("Done.")
def replaceAccession(out):

    with open(out,'r') as infile:
        locus = ""
        for line in infile:

            if line.startswith("LOCUS"):
                locus = line.split()[1]
                print(line.rstrip())
            elif line.startswith("ACCESSION"):
                l = line.split()
                l[1] = locus
                #print("ACCESSION {0:<3}".format(locus))
                print('   '.join(l))
            elif line.startswith("VERSION"):
                l = line.split()
                l[1] = locus
                print('     '.join(l) )
            else:
                print(line.rstrip())


def main():
    parser = argparse.ArgumentParser(description="Script to replace locus [contig name] tag in genbank records using a tab delimited mapping file (may be generated with rename fasta headsers script).")
    parser.add_argument(dest='genbank_file', metavar='genbankFile', help='Genbank file to process.')
    parser.add_argument(dest='auxiliar_file', metavar='auxiliarFile', help='2 columns tab delimited file matching old locus ID [col 0] to the new ones to replace [col 1].')
    parser.add_argument(dest='output_file', metavar='outputFile', help='Output file.')
    args = parser.parse_args()

    dict = processAuxFile(args.auxiliar_file)
    processGenbank(args.genbank_file,args.output_file,dict)
    replaceAccession(args.output_file)
if __name__ == '__main__':
    main()
