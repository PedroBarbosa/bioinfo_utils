from Bio import SeqIO
import argparse
import sys


def readInput(infile,outbasename,window):
    bedout=outbasename + "_exactAdaptorCoordinates.bed"
    bedoutWindow=outbasename + "_extended" + str(window) + "CoordinatesAroundAdaptors.bed"
    with open(infile,'r') as infl:
        with open(bedout,'w') as outbed:
            with open(bedoutWindow,'w') as outbedextended:
                trim=False
                for line in infl:

                    if trim and not "Sequence name" in line and line.strip():
                        fields=line.rstrip().split("\t")
                        if ',' in fields[2]:
                            coords = fields[2].split(",")
                            for i in coords:
                                if int(i.split("..")[0]) == 1 or int(i.split("..")[1]) == int(fields[1]):
                                    outbed.write("%s\t%s\t%s\n" % (fields[0],int(i.split("..")[0])-1,i.split("..")[1]))
                                    outbedextended.write("##%s\t%s\t%s\t%s\t%s\n" % (fields[0],fields[1],i,fields[3],"TERMINAL"))
                                else:
                                    outbed.write("%s\t%s\t%s\n" % (fields[0],int(i.split("..")[0])-1,i.split("..")[1]))
                                    outbedextended.write("%s\t%i\t%i\n" % (fields[0],int(i.split("..")[0]) - window -1, int(i.split('..')[1]) + window))
                        else:
                            if int(fields[2].split("..")[0]) == 1 or int(fields[2].split("..")[1]) == int(fields[1]):
                                outbed.write("%s\t%s\t%s\n" % (fields[0],int(fields[2].split("..")[0])-1,fields[2].split("..")[1]))
                                outbedextended.write("##" + "\t".join(fields) + "\t" + "TERMINAL\n")
                            else:
                                outbed.write("%s\t%s\t%s\n" % (fields[0],int(fields[2].split("..")[0]) -1,fields[2].split("..")[1]))
                                outbedextended.write("%s\t%i\t%i\n" % (fields[0], int(fields[2].split("..")[0]) - window -1, int(fields[2].split("..")[1]) + window))
                    if line.startswith("Trim"):
                        trim=True
def main():
    parser = argparse.ArgumentParser(description='Script to generate bed file to inspect more carefully the regions flagged to contain sequence adpators using as inout the contamination file produced by NCBI screening pipeline.')
    parser.add_argument(dest='inputfile', metavar='contaminationReport', help='Input file.')
    parser.add_argument(dest='outbasename', metavar='outbasename', help='Output basename')
    parser.add_argument("-w", metavar='window', type=int, default=20,help="Number of base pairs to extend genome intervals flagged as adaptors. Default: 20")
    args = parser.parse_args()

    readInput(args.inputfile,args.outbasename,args.w)
if __name__ == "__main__":
    main()