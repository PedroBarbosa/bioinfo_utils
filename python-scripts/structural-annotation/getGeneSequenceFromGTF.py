__author__ = 'pedro'

import argparse
import gffutils
from gffutils import helpers
from collections import defaultdict
import os
import logging
import sys
from Bio import SeqIO
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')


def extractFromFasta(reference_fasta, dict_plus, dict_minus, dict_unknown, output_basename):
    logging.info("Processing reference fasta file and writing output...")
    if os.path.exists(output_basename + '.fasta'):
        os.remove(output_basename + '.fasta')

    if os.path.exists(output_basename + '_unknownStrand.fasta'):
        os.remove(output_basename + '_unknownStrand.fasta')

    if len(dict_unknown) == 0:
        logging.info("Good!! No genes with unknown strand information were found. All the genes will be written to the output file with proper orientation.")

    else:
        logging.info("%i genes have unknown strand information. They will be written to a different file based on coordinates of the scaffold they belong." % len(dict_unknown))

    fasta_records = 0
    scaffolds_with_no_genes = 0

    with open(output_basename + '.fasta', "w") as file:
        with open(output_basename + '_unknownStrand.fasta', "w") as file_unknown:
            handle = open(reference_fasta, "rU")
            for record in SeqIO.parse(handle, "fasta"):
                fasta_records += 1
                if record.id in dict_plus:
                    genes_per_scaffold = dict_plus[record.id]
                    for gene in genes_per_scaffold:
                        file.write(">" + gene[0] + "\n")
                        file.write(str(record.seq[int(gene[1] - 1): int(gene[2] - 1)]) + "\n")

                if record.id in dict_minus:
                    genes_per_scaffold = dict_minus[record.id]
                    for gene in genes_per_scaffold:
                        file.write(">" + gene[0] + "\n")
                        sequence = record.seq[int(gene[1] - 1): int(gene[2] - 1)]
                        file.write(str(sequence[::-1]) + "\n") #reverse sequence

                if len(dict_unknown) > 0 and record.id in dict_unknown:
                    genes_per_scaffold = dict_unknown[record.id]
                    for gene in genes_per_scaffold:
                        file_unknown.write(">" + gene[0] + "\n")
                        file_unknown.write(str(record.seq[int(gene[1] - 1): int(gene[2] - 1)]) + "\n")


                if not record.id in dict_plus and not record.id in dict_minus and not record.id in dict_unknown:
                    scaffolds_with_no_genes += 1

            if os.stat(output_basename + '_unknownStrand.fasta').st_size == 0:
                os.remove(output_basename + '_unknownStrand.fasta')

    logging.info("Info: %i scaffolds out of %i have no annotated genes in the GTF file." % (scaffolds_with_no_genes,fasta_records))
    file.close()
    logging.info("Done.")



def createDicGenesPerScaffold(db):
    logging.info("Creating dictionary of genes per scaffold..")
    dict_plus = defaultdict(list)
    dict_minus = defaultdict(list)
    dict_unknown = defaultdict(list)

    for gene in db.features_of_type('gene',strand='+') :
        scaffold_id = gene.seqid
        dict_plus[scaffold_id].append((gene.id,gene.start,gene.end))

    for gene in db.features_of_type('gene',strand='-') :
        scaffold_id = gene.seqid
        dict_minus[scaffold_id].append((gene.id,gene.start,gene.end))

    for gene in db.features_of_type('gene',strand='.') :
        scaffold_id = gene.seqid
        dict_unknown[scaffold_id].append((gene.id,gene.start,gene.end))


    return dict_plus,dict_minus, dict_unknown


def createGffUtilsGffread(gtf_file):
    dbname=os.path.basename(gtf_file).split('.')[0] + "DB.sql"
    try:
        db=gffutils.FeatureDB(dbname)
        logging.info("Previously generated database successfully loaded")
    except:
        logging.info("Creating new gffutils database..")
        dialect=helpers.infer_dialect(['PE100bp.genome1_contig-1000002  AUGUSTUS        CDS     7072    7639    0.47    +       0       transcript_id "PE100bp.genome1_contig-1000002.g78050.t1"; gene_id "PE100bp.genome1_contig-1000002.g78050";'])
        db=gffutils.create_db(gtf_file, dbfn=dbname, id_spec={'gene': ['gene_id', 'gene_name', 'nearest_ref'],'transcript' : ['transcript_id', 'oId'], 'exon': 'exon_number'},
                          merge_strategy="create_unique",keep_order=True, sort_attribute_values=True, disable_infer_transcripts=False, disable_infer_genes=False,
                          dialect=dialect, checklines=50 ,verbose=False,force=True)
        logging.info("Database " + dbname + " sucessfully generated.")
        return db

def createGffUtilsCuffmerge(gtf_file):
    dbname=os.path.basename(gtf_file).split('.')[0] + "DB.sql"
    try:
        db=gffutils.FeatureDB(dbname)
        logging.info("Previously generated database successfully loaded.")
    except:
        logging.info("Creating new gffutils database...")
        dialect=helpers.infer_dialect(['Potrx000002 Cufflinks       exon    8052    8625    .       -       .       gene_id "XLOC_000003"; transcript_id "TCONS_00000005"; exon_number "2"; '
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
    parser.add_argument(dest='output_basename', metavar='output_basename', nargs=1, help='Basename to which the output files will be written.')
    parser.add_argument(dest='software', metavar='software', nargs=1, type=str, choices=['cuffmerge','gffread'],help='Tool that produced the input gtf file.')

    args = parser.parse_args()

    if "cuffmerge" in args.software:
        db = createGffUtilsCuffmerge(args.gtf_file[0])
    elif "gffread" in args.software:
        db = createGffUtilsGffread(args.gtf_file[0])
    dicScaffGenes_plus, dicScaffGenes_minus, dicScaffGenes_unknown =  createDicGenesPerScaffold(db)
    extractFromFasta(args.reference_fasta[0], dicScaffGenes_plus, dicScaffGenes_minus, dicScaffGenes_unknown, args.output_basename[0])
    
     
if __name__ == "__main__":
    main()
