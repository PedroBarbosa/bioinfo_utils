import argparse
from Bio import SeqIO
from collections import defaultdict

def replaceFasta(infasta,inbed,outfasta):
    with open(inbed,'r') as inb:
        dict_bed= defaultdict(list)
        for line in inb:
            if not line.startswith("#"):
                fields=line.rstrip().split("\t")
                dict_bed[fields[0]].append((fields[1],fields[2]))
    inb.close()
    for k, v in dict_bed.items():
        print(k,v)


    handle = open(infasta,"rU")
    handleout = open(outfasta,"w")
    for seq_record in SeqIO.parse(handle,"fasta"):
        header = seq_record.description
        if header not in dict_bed.keys():
            handleout.write(header + "\n" + seq_record.seq + "\n")
        else:
            seqToEdit=seq_record.seq
            for adaptor in dict_bed[header]:
                seqToEdit = seqToEdit.replace(seqToEdit[adaptor[0]:adaptor[1]],"N")
                print(seqToEdit)


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
