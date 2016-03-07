import argparse
import operator

def processFastaFiles(inputFile):

    contig_seq = ""
    contig_id = ""
    final_dict = {}
    with open(inputFile[0]) as file:
        for line in file:
            line=line.rstrip()
            if line.startswith('>'):

                if contig_id:

                    previous_id = contig_id
                    final_dict[previous_id] = len(contig_seq)

                contig_id = line
                contig_seq = ""


            else:
                contig_seq += line



        ###process last contig (after last >) ###
        final_dict[contig_id] = len(contig_seq)

    file.close()
    return  final_dict


def writeOutput(final_dict, outputFile):
    sorted_dict = sorted(final_dict.items(), key=operator.itemgetter(1), reverse=True)
    with open(outputFile, "w") as out:
        for k,v in sorted_dict.iteritems():
            out.write(k + " " + v + "\n")


parser = argparse.ArgumentParser(description='Script to create genome table from fasta file. [chromossome_id length_bp].')
parser.add_argument(dest='inputFile', metavar='input_files', nargs=1,
                    help='Fasta file to be processed.')
parser.add_argument('-o', '--outputFile', metavar='', required = True, help='Output file')
args = parser.parse_args()


dict = processFastaFiles(args.input_files[0])
writeOutput(dict,args.outputFile)