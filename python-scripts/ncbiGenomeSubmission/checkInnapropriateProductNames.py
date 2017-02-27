import argparse
from Bio import SeqIO


def processGenbank(gbk, out, listlocus):
    list = [line.rstrip() for line in open(listlocus)]
    with open(out,'w') as outfile:
        for gb_record in SeqIO.parse(open(gbk, "r"), "genbank"):
            for feature in gb_record.features:
                if feature.type == "CDS" and feature.qualifiers["locus_tag"][0] in list:
                    if "/M" in feature.qualifiers['product'][0]:
                        new_product=feature.qualifiers['product'][0].split("/")[0]
                        feature.qualifiers['product'] = new_product

            SeqIO.write(gb_record, outfile, "genbank")

def main():
    parser = argparse.ArgumentParser(
        description="Script to replace locus [contig name] tag in genbank records using a tab delimited mapping file (may be generated with rename fasta headsers script).")
    parser.add_argument(dest='genbank_file', metavar='genbankFile', help='Genbank file to process.')
    parser.add_argument(dest='locus_tag_list', metavar='locusTagList', help='File listing locus tag to have a deeper look. One per line.')
    parser.add_argument(dest='output_file', metavar='outputFile', help='Output file.')
    args = parser.parse_args()

    processGenbank(args.genbank_file, args.output_file, args.locus_tag_list)


if __name__ == '__main__':
    main()
