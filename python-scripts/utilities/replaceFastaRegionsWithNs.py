import argparse
from Bio import SeqIO
from collections import defaultdict

def replaceFasta(infasta,inbed,outfasta):
    with open(inbed,'r') as inb:
        dict_bed= defaultdict(list)
        for line in inb:
            if not line.startswith("#"):
                fields=line.rstrip().split("\t")
                dict_bed[fields[0]].append((int(fields[1]),int(fields[2])))
    inb.close()
    handle = open(infasta,"rU")
    handleout = open(outfasta,"w")
    before=0
    after=0
    noNs=0
    for seq_record in SeqIO.parse(handle,"fasta"):
        header = seq_record.description
        ed=""
        if header not in dict_bed.keys():
            handleout.write(">" + header + "\n" + str(seq_record.seq) + "\n")
        else:
            i=1
            for adaptor in dict_bed[header]:
                Ns = "N" * (adaptor[1]- adaptor[0])
                edited = "".join((str(seq_record.seq[:adaptor[0]]),Ns,str(seq_record.seq[adaptor[1]:])))
                upstream = int(adaptor[0] - 5)
                downstream = int(adaptor[1] + 5)
                if "NNNNN" in seq_record.seq[upstream:adaptor[0]]:
                    before+=1
                elif "NNNNN" in seq_record.seq[adaptor[1]:downstream]:
                    after+=1
                else:
                    noNs+=1
                print("".join((header, "     ",str(seq_record.seq[upstream:adaptor[0]]), "     ",Ns, "     ", str(seq_record.seq[adaptor[1]:downstream]))))

                #Process sequences for which there are multiple bed regions to replace
                if len(dict_bed[header]) > 1 and i==1:
                    i+=1
                    ed="".join((edited[:adaptor[0]],Ns,edited[adaptor[1]:]))
                elif len(dict_bed[header]) > 1 and len(dict_bed[header]) == i:
                    ed="".join((ed[:adaptor[0]],Ns,ed[adaptor[1]:]))
                    handleout.write(">" + header + "\n" + ed + "\n")
                elif len(dict_bed[header]) > 1 and i < len(dict_bed[header]):
                    i+=1
                    ed = "".join((ed[:adaptor[0]],Ns,ed[adaptor[1]:]))
                elif len(dict_bed[header]) == 1:
                    handleout.write(">" + header + "\n" + edited + "\n")

    print("%s:%i\n%s:%i\n%s:%i" % ("Regions with 5 continuous Ns before the region replaced",before,"Regions with 5 continuous Ns after the region replaced",
                                   after, "Regions with no 5 Ns surrounding the replaced region",noNs))

def main():
    parser = argparse.ArgumentParser(description='Script to replace stretches of a genome with undetermined nucleotides (Ns) based on a bed file. [Done this script'
                                                 'to correct contaminations with adaptors in the ncbi submission process of the cork oak draft geneome]')
    parser.add_argument(dest='fastafile', metavar='fastaFile', help='Input fasta file.')
    parser.add_argument(dest='bedfile', metavar='bedFile', help='Input bed file.')
    parser.add_argument(dest='outfile', metavar='outFile', help='Output file.')
    args = parser.parse_args()
    replaceFasta(args.fastafile,args.bedfile,args.outfile)

if __name__ == "__main__":
    main()

