import argparse
from Bio import SeqIO


def processGenbank(gbk, out, prefix):
    i = "00001"
    print("Replacing genbank records..")
    with open(out, 'w') as outfile:
        for gb_record in SeqIO.parse(open(gbk, "r"), "genbank"):
            for feature in gb_record.features:
                if feature.type == "gene":
                    feature.qualifiers["locus_tag"] = prefix + "_" + i
                    i=str(int(i) + 1).zfill(len(i))

                elif feature.type != "source":
                    feature.qualifiers["locus_tag"] = prefix + "_" + str(int(i) - 1).zfill(len(i))


            SeqIO.write(gb_record, outfile, "genbank")

def main():
    parser = argparse.ArgumentParser(
        description="Script to replace locus [contig name] tag in genbank records using a tab delimited mapping file (may be generated with rename fasta headsers script).")
    parser.add_argument(dest='genbank_file', metavar='genbankFile', help='Genbank file to process.')
    parser.add_argument(dest='locus_tag_prefix', metavar='locusTagPrefix', help='String to represent locus tag prefix.')
    parser.add_argument(dest='output_file', metavar='outputFile', help='Output file.')
    args = parser.parse_args()

    processGenbank(args.genbank_file, args.output_file, args.locus_tag_prefix)


if __name__ == '__main__':
    main()
