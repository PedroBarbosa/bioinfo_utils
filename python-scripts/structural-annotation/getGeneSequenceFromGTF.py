__author__ = 'pedro'

import argparse
import gffutils
from gffutils import helpers
import os
import logging
import sys
from Bio import SeqIO
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')


def extractFromFasta(reference_fasta, dict, output_file):
    logging.info("Processing reference fasta file and writing output...")
    if os.path.exists(output_file):
        os.remove(output_file)

    scaffolds_with_no_genes = 0
    with open(output_file, "w") as file:
        handle = open(reference_fasta, "rU")
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in dict:
                genes_per_scaffold = dict[record.id]
                for gene in genes_per_scaffold:
                    file.write(">" + gene[0] + "\n")
                    file.write(str(record.seq[int(gene[1] - 1): int(gene[2] - 1)]) + "\n")

            else:
                scaffolds_with_no_genes += 1

    logging.info("Info:" + str(scaffolds_with_no_genes) + " scaffolds have no annotated genes in the GTF file.")
    file.close()
    logging.info("Done.")



def createDicGenesPerScaffold(db):
    logging.info("Creating dictionary of genes per scaffold..")
    dict = {}
    for gene in db.features_of_type('gene'):
        gene_id = gene.id
        start = gene.start
        end = gene.end
        gene_tuple = (gene_id,start,end)
        scaffold_id = gene.seqid
        if scaffold_id in dict:
            new_list = dict[scaffold_id]
            new_list.append(gene_tuple)
            dict[scaffold_id] = new_list
        else:
            dict[scaffold_id] = [gene_tuple]

    return dict

def createGffUtilsCuffmerge(gtf_file):

    logging.info("Creating gffutils database...")
    dbname=os.path.abspath(gtf_file).split('.')[0] + "DB.sql"

    dialect=helpers.infer_dialect(['Potrx000002	Cufflinks	exon	8052	8625	.	-	.	gene_id "XLOC_000003"; transcript_id "TCONS_00000005"; exon_number "2"; '
                                   'gene_name "Potrx000002g00030"; oId "Potrx000002g00030.1"; nearest_ref "Potrx000002g00030.1"; class_code "="; tss_id "TSS3"; p_id "P3";'])

    db=gffutils.create_db(gtf_file, dbfn=dbname, id_spec={'gene': ['gene_id', 'gene_name', 'nearest_ref'],'transcript' : ['transcript_id', 'oId'], 'exon': 'exon_number'},
                          merge_strategy="create_unique",keep_order=True, sort_attribute_values=True, disable_infer_transcripts=False, disable_infer_genes=False,
                          dialect=dialect, checklines=50 ,verbose=False,force=True)
    logging.info("Database " + dbname + " successfuly generated.")

    return db


def main():
    parser = argparse.ArgumentParser(description='Script to extract the gene sequences from a genome reference file and its corresponding GTF. Requires gffutils and Biopython modules to be installed.')
    parser.add_argument(dest='gtf_file', metavar='gtf_file', nargs=1, help='GTF file to be processed.')
    parser.add_argument(dest='reference_fasta', metavar='reference_fasta', nargs=1, type=str,help='FASTA file representing the genome sequence.')
    parser.add_argument(dest='output_file', metavar='output_file', nargs=1, help='File to write the output.')
    parser.add_argument(dest='software', metavar='software', nargs=1, type=str, choices=['cuffmerge'],help='Tool that produced the input gtf file.')

    args = parser.parse_args()

    if "cuffmerge" in args.software:

        db = createGffUtilsCuffmerge(args.gtf_file[0])
        dicScaffGenes =  createDicGenesPerScaffold(db)
        extractFromFasta(args.reference_fasta[0], dicScaffGenes, args.output_file[0])

if __name__ == "__main__":
    main()