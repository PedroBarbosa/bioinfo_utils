__author__ = 'pedro'

import argparse
import os
import gffutils




def createGffUtilsCuffmerge(gtf_file):
    dbname=os.path.basename(gtf_file).split('.')[0]
    gffutils.create_db(gtf_file, dbfn=dbname, id_spec={'gene': ['gene_id', 'gene_name'],'transcript' : ['transcript_id']},)

def main():
    parser = argparse.ArgumentParser(description='Script to check the relationship between genes, transcripts and exons in a GTF file.')
    parser.add_argument(dest='gtf_file', metavar='gtf_file', nargs=1, help='GTF file to be processed.')
    parser.add_argument(dest='software', metavar='software', nargs=1, type=str, choices=['cuffmerge','cufflinks','stringtie'],help='Tool that produced the input gtf file.')
    #parser.add_argument(dest='output_prefix', metavar='output_prefix', nargs=1, help='Basename to write on the ouptut files.')
    args = parser.parse_args()


    if "cuffmerge" in args.software:
        createGffUtilsCuffmerge(args.gtf_file[0])

if __name__ == "__main__":
    main()