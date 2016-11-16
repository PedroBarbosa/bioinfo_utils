import argparse
import os

def processFastaFiles(inputFile):

    contig_seq = ""
    contig_id = ""
    final_dict = {}
    ordered_fasta = {}
    with open(inputFile) as file:
        for line in file:
            line=line.rstrip()
            if line.startswith('>'):

                if contig_id:

                    previous_id = contig_id
                    final_dict[previous_id[1:]] = (len(contig_seq),contig_seq)
                contig_id = line
                contig_seq = ""

            else:
                contig_seq += line

        ###process last contig (after last >) ###
        final_dict[contig_id[1:]] = (len(contig_seq),contig_seq)

    file.close()
    return  final_dict


def writeOutput(final_dict,orderedFasta,outputFile, inputFileName):
    sorted_dict = sorted(final_dict.keys(), key=lambda x: final_dict[x], reverse=True)
    with open(outputFile, "w") as out:
        for k,v in iter(sorted_dict.items()):
            out.write(k + "\t" + str(v[0]) + "\n")

    if orderedFasta:
        outFastaFile = os.getcwd() +'/' + str(os.path.basename(inputFileName).rsplit('.')[0]) + "-ordered.fasta"
        with open(outFastaFile, "w") as outFasta:
            for k,v in iter(sorted_dict.items()):
                outFasta.write(k + "\n" + str(v[1]) + "\n")


parser = argparse.ArgumentParser(description='Script to create genome table from fasta file. [chromossome_id length_bp].')
parser.add_argument(dest='inputFile', metavar='input_files', nargs=1,
                    help='Fasta file to be processed.')
parser.add_argument('-o', '--outputFile', metavar='', required = True, help='Output file')
parser.add_argument('-f' '--fastaOrdered',action='store_true', help='Flag that writes new fasta ordered by legnth of sequences.' )
args = parser.parse_args()


dict_out = processFastaFiles(args.inputFile[0])
writeOutput(dict_out,args.fastaOrdered,args.outputFile,args.inputFile[0])