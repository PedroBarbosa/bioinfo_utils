__author__ = 'pedro'

import argparse
import os
import gffutils
from gffutils import helpers
import logging

global number_of_genes
global number_of_transcripts
global number_of_exons

def generalStats(db):
    number_of_genes = db.count_features_of_type("gene")
    number_of_transcripts=db.count_features_of_type("transcript")
    number_of_exons=db.count_features_of_type("exon")

def transcriptExonUsage(db):


    exons= db.features_of_type("exon")

    for gene in db.features_of_type("gene"):

        transcripts = db.children(gene.id, level=1, featuretype='transcript')
        for transcript in transcripts:
            #write bed to file and perform merge within python with subprocess module
            print db.bed12(transcript)


        #generate several stats according the tutorial
        merged_transcripts=db.merge(transcripts)
        for merged_transcript in merged_transcripts:
            print merged_transcript.attributes
     #   fout = open('/home/pedro/Desktop/new.gtf','w')
      #  for merged_exon in merged_exons:
       #     print [transcript for transcript in db.parents(merged_exon,1,featuretype='transcript')]


            #[transcript.id for transcript.id in db.children(gene.id, featuretype='transcript')]
        #for transcript in db.children(gene.id, featuretype="transcript"):
            #print gene.id, transcript.id
       # [list_exons.id for list_exons in db.children(transcript.id, featuretype="exon")]

        print "\n\n\n"



def createGffUtilsCuffmerge(gtf_file):
    dbname=os.path.basename(gtf_file).split('.')[0] + "DB"
    dialect=helpers.infer_dialect(['gene_id "XLOC_000001"; transcript_id "TCONS_00000001"; exon_number "1"; gene_name "Potrx000002g00010"; oId "Potrx000002g00010.1"; '
                                   'nearest_ref "Potrx000002g00010.1"; class_code "="; tss_id "TSS1"; p_id "P1";'])

    try:
        db=gffutils.create_db(gtf_file, dbfn=dbname, id_spec={'gene': ['gene_id', 'gene_name', 'nearest_ref'],'transcript' : ['transcript_id', 'oId'], 'exon': 'exon_id'},merge_strategy="merge",
                              disable_infer_transcripts=False, disable_infer_genes=False, dialect=dialect, checklines=1000 ,verbose=True,force=True)
    except:
        logging.warning("Database already exists. No need to create new one.")
        db=gffutils.FeatureDB(dbname)


    generalStats(db)
    transcriptExonUsage(db)

    #     print gene_id.id
    #for gene in db.children(gene_id.id,order_by='featuretype'):#="exon"):
    #         print'{0.featuretype:>12}: {0.id}'.format(gene)

  #  for gene in db.featuretypes():
 #       print db.children(gene)
    #       print i.bin
    #      print i.dialect
#print "\n\n"

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