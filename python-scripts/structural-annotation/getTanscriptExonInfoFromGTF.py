__author__ = 'pedro'

import argparse
import os
import gffutils
from gffutils import helpers
import logging
import re
import subprocess


global number_of_genes
global number_of_transcripts
global number_of_exons

def checkBedtoolsInstallation(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def generalStats(db):
    number_of_genes = db.count_features_of_type("gene")
    number_of_transcripts=db.count_features_of_type("transcript")
    number_of_exons=db.count_features_of_type("exon")

def executeGeneTranscripExonUsage(exons_string,exonsPertranscript,gene2transcript,gene_id,transcript_id,exon_id, repeated_transcripts):

#    if gene_id not in gene2transcript and len(gene2transcript)!=0: #Process last transcript when new gene appears avoiding also the first line of file
        #Process the last transcript
#        if exons_string not in exonsPertranscript: #While set of exons per transcript is new, add to the list of exons set
#            exonsPertranscript.append(exons_string)

#        else: #If set of exon per transcript is repeated, transcript is repeated
#            repeated_transcripts += 1





    if transcript_id in gene2transcript[gene_id]: #While in the same transcript, populate the exons string
        exons_string += exon_id + ";"

    else: #New transcript
        #Process the last transcript
        if exons_string not in exonsPertranscript: #While set of exons per transcript is new, add to the list of exons set
            exonsPertranscript.append(exons_string)

        else: #If set of exon per transcript is repeated, transcript is repeated
            repeated_transcripts += 1

            #Process the actual one
        gene2transcript[gene_id].append(transcript_id) #update list of transcripts per gene
        exons_string = exon_id + ";"


    return exons_string,exonsPertranscript,gene2transcript,repeated_transcripts

def transcriptWithSameExonsCuffmerge(gtffile):
    gene2transcript= {}
    exonsPertranscript = []
    repeated_transcripts = 0
    exons_string=""
    with open(gtffile,"r") as file:
        for line in file:

            feature = re.split(r'\t',line)
            attributes = re.split(r';',feature[8])
            gene_id = re.split(r'\"',attributes[0])[1]
            transcript_id = re.split(r'\"',attributes[1])[1]
            exon_id = re.split(r'\"',attributes[2])[1]


            if gene_id in gene2transcript:
                exons_string,exonsPertranscript,gene2transcript,repeated_transcripts = executeGeneTranscripExonUsage(exons_string,exonsPertranscript,gene2transcript,gene_id,transcript_id,exon_id,repeated_transcripts)


#            elif len(gene2transcript) != 0:
            else:
                gene2transcript[gene_id] = [transcript_id] #create new gene transcript instance
                if exons_string not in exonsPertranscript: #While set of exons per transcript is new, add to the list of exons set
                    exonsPertranscript.append(exons_string)

                else: #If set of exon per transcript is repeated, transcript is repeated
                    repeated_transcripts += 1
                    
                exons_string,exonsPertranscript,gene2transcript,repeated_transcripts = executeGeneTranscripExonUsage(exons_string,exonsPertranscript,gene2transcript,gene_id,transcript_id,exon_id,repeated_transcripts)

                #exonsPertranscript = []
                #exons_string=""
                #exons_string,exonsPertranscript,gene2transcript,repeated_transcripts = executeGeneTranscripExonUsage(exons_string,exonsPertranscript,gene2transcript,gene_id,transcript_id,exon_id,repeated_transcripts)



            # else: #First line
            #     gene2transcript[gene_id] = [transcript_id]
            #     exons_string = exon_id + ";"





def transcriptExonUsage(db):

    fout = open('/home/pedro/Desktop/new.bed','w')
    exons= db.features_of_type("exon")

    for gene in db.features_of_type("gene"):

        transcripts = db.children(gene.id, level=1, featuretype='transcript')
        for transcript in transcripts:
            #write bed to file and perform merge within python with subprocess module
            fout.write(db.bed12(transcript)+"\n")
   #         ps = subprocess.Popen(['grep', '>', fastaReadsFile], stdout=subprocess.PIPE)
    #return subprocess.check_output(['wc', '-l'], stdin=ps.stdout)
    #ps.wait()
#    print checkBedtoolsInstallation("bedtools")

        #generate several stats according the tutorial
  #      merged_transcripts=db.merge(transcripts)
  #      for merged_transcript in merged_transcripts:
  #          print merged_transcript.attributes
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
                              disable_infer_transcripts=False, disable_infer_genes=False, dialect=dialect, checklines=1000 ,verbose=True,force=False)
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
        #createGffUtilsCuffmerge(args.gtf_file[0])
        transcriptWithSameExonsCuffmerge(args.gtf_file[0])
if __name__ == "__main__":
    main()