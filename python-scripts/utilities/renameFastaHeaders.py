from Bio import SeqIO
import argparse
import sys

parser = argparse.ArgumentParser(description='Script to rename the header of a FASTA file, either by adding string in the end of header or totally replace it. Biopython needed')
parser.add_argument(dest='fastaFile', metavar='fastaFile', help='FASTA file to process.')
parser.add_argument('-a', metavar='add', help='String to add to the end of each header. String will be concatenated to each header with no spaces.')
parser.add_argument('-r', metavar='replace', help='Base string to replace all headers completely. For each record a "_i" number will be incremented.')
parser.add_argument('-w', '--removeSpace', action='store_true',help='Flag to select first string up to the first whitespace.')
parser.add_argument('-s', metavar='split', help='Split headers by the last ocurrence of given char and select leading substring. If char to split not found, old header will be kept.')
parser.add_argument('-f', metavar='subsetFile', help='Only rename and output headers present in this file [One per line].' )
args = parser.parse_args()

subsetList=[]
if args.f:
    with open(args.f, 'r') as infile:
        for line in infile:
            subsetList.append(line.rstrip())

if args.a and args.r and args.s and args.removeSpace:
    raise ValueError("Choose only one of the options: add [-a], replace [-r], split [-s] or  remove [-w]")

elif args.a:

    handle = open(args.fastaFile,"rU")
    handle2 = open("correspondingSequenceNames.tsv","w")

    for seq_record in SeqIO.parse(handle,"fasta"):
        old_header = seq_record.description
        seq_record.id = old_header + args.a
        if len(subsetList) > 0:
            if old_header in subsetList:
                print('>' + seq_record.id + '\n' + seq_record.seq)
                handle2.write(old_header+ '\t' + seq_record.id +'\n')
        else:
            print('>' + seq_record.id + '\n' + seq_record.seq)
            handle2.write(old_header+ '\t' + seq_record.id +'\n')

    handle.close()
    handle2.close()

elif args.r:
    i=1
    handle = open(args.fastaFile,"rU")
    handle2 = open("correspondingSequenceNames.tsv","w")
    for seq_record in SeqIO.parse(handle,"fasta"):
        old_header = seq_record.description
        seq_record.id = args.r + '_' + str(i)
        if len(subsetList) > 0:
            if old_header in subsetList:
                print('>' + seq_record.id + '\n' + seq_record.seq)
                handle2.write(old_header+ '\t' + seq_record.id +'\n')
                i+=1
        else:
            print('>' + seq_record.id + '\n' + seq_record.seq)
            handle2.write(old_header+ '\t' + seq_record.id +'\n')
            i+=1
    handle.close()
    handle2.close()

elif args.s:
    handle=open(args.fastaFile,"rU")
    handle2 = open("correspondingSequenceNames.tsv","w")
    for seq_record in SeqIO.parse(handle,"fasta"):
        old_header = seq_record.description
        seq_record.id = old_header.rsplit(args.s,1)[0]
        if len(subsetList) > 0:
            if old_header in subsetList:
                if seq_record.id:
                     print('>' + seq_record.id + '\n' + seq_record.seq)
                     handle2.write(old_header+ '\t' + seq_record.id +'\n')
                else:
                     print('>' + old_header + '\n' + seq_record.seq)
                     handle2.write(old_header+ '\t' + old_header +'\n')
        else:
            if seq_record.id:
                print('>' + seq_record.id + '\n' + seq_record.seq)
                handle2.write(old_header+ '\t' + seq_record.id +'\n')
            else:
                print('>' + old_header + '\n' + seq_record.seq)
                handle2.write(old_header+ '\t' + old_header +'\n')
    handle.close()

elif args.removeSpace:
    handle=open(args.fastaFile,"rU")
    handle2 = open("correspondingSequenceNames.tsv","w")
    for seq_record in SeqIO.parse(handle,"fasta"):
        old_header = seq_record.description
        seq_record.id = old_header.split(" ")[0]
        if len(subsetList) > 0:
            if old_header in subsetList:
                print('>' + seq_record.id + '\n' + seq_record.seq)
                handle2.write(old_header+ '\t' + seq_record.id +'\n')
        else:
            print('>' + seq_record.id + '\n' + seq_record.seq)
            handle2.write(old_header+ '\t' + seq_record.id +'\n')
    handle.close()

else:
    raise ValueError("Please set one of the options: add [-a] or replace [-r] or splot [-s] or remove [-w]")
