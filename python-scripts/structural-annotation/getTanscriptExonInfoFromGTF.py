__author__ = 'pedro'

import argparse
import os
import gffutils
from gffutils import helpers
import logging
import re
import sys
import csv


global number_of_genes
global number_of_transcripts
global number_of_exons








def executeGeneTranscripExonUsage(exons_string,exonsPertranscript,gene2transcript,gene_id,transcript_id,exon_id, numb_repeated_transcripts,previous_transcript_id, repeated_transcritps):

    if gene_id not in gene2transcript: #Process last transcript when new gene appears

        if exons_string in exonsPertranscript: #Since is last transcript of gene, no need to add to the set of exons per transcript
            numb_repeated_transcripts += 1
            repeated_transcritps.append(previous_transcript_id)
        gene2transcript[gene_id] = [transcript_id]
        exons_string = ""
        exonsPertranscript = []


    if transcript_id in gene2transcript[gene_id]: #While in the same transcript, populate the exons string
        exons_string += exon_id + ";"

    else: #New transcript
        #Process the last transcript
        if exons_string not in exonsPertranscript: #While set of exons per transcript is new, add to the list of exons set
            exonsPertranscript.append(exons_string)

        else: #If set of exon per transcript is repeated, transcript is repeated
            numb_repeated_transcripts += 1
            repeated_transcritps.append(previous_transcript_id)

        #Process the actual one
        gene2transcript[gene_id].append(transcript_id) #update list of transcripts per gene
        exons_string = exon_id + ";"

    return exons_string,exonsPertranscript,gene2transcript,numb_repeated_transcripts,repeated_transcritps



def getRepeatedTranscripts(gtffile):
    gene2transcript= {}
    exonsPertranscript = []
    numb_repeated_transcripts = 0
    repeated_transcripts = []
    exons_string=""
    transcript_id = ""
    original_numb_exon= 0
    with open(gtffile,"r") as file:
        for line in file:
            if not line.startswith("#"):
                line.rstrip()
                feature = re.split(r'\t',line)
                try:
                    if feature[2] == "exon":
                        original_numb_exon += 1
                        attributes = re.split(r';',feature[8])
                        gene_id = re.split(r'\"',attributes[0])[1]
                        previous_transcript_id = transcript_id
                        transcript_id = re.split(r'\"',attributes[1])[1]
                        exon_id = re.split(r'\"',attributes[2])[1]


                        if gene_id in gene2transcript:
                            exons_string,exonsPertranscript,gene2transcript,numb_repeated_transcripts,repeated_transcripts = executeGeneTranscripExonUsage(exons_string,exonsPertranscript,gene2transcript,gene_id,transcript_id,exon_id,numb_repeated_transcripts,previous_transcript_id, repeated_transcripts)

                        elif len(gene2transcript) != 0:
                            exons_string,exonsPertranscript,gene2transcript,numb_repeated_transcripts,repeated_transcripts = executeGeneTranscripExonUsage(exons_string,exonsPertranscript,gene2transcript,gene_id,transcript_id,exon_id,numb_repeated_transcripts,previous_transcript_id, repeated_transcripts)

                        else: #First line
                            gene2transcript[gene_id] = [transcript_id]
                            exons_string =  exon_id + ";"
                except IndexError:
                    sys.stderr.write("ERROR with input file.Please check your gtf file, might exist some empty lines at the end of the file.\n")
                    sys.exit(2)

        if exons_string in exonsPertranscript: #Last line of file
            numb_repeated_transcripts += 1
            repeated_transcripts.append(transcript_id) #not previous because file is read


    return numb_repeated_transcripts, repeated_transcripts,original_numb_exon



def newGTFmerged(db,outputbasename,software):
    if os.path.exists(outputbasename + 'merged.gtf'):
        os.remove(outputbasename + 'merged.gtf')

    def newattribute(feature,newAttributes):
        newAttributes.pop()
        new_attributes_string=""
        for att in newAttributes:
            new_attributes_string += att + ";"

        feature[8] = new_attributes_string #assign new attributes string
        return '\t'.join(feature) + '\n'

    def gene_id_rearrangement_in_transcripts(GffutilsTranscriptFeature):
        transcript_feature = re.split(r'\t',str(GffutilsTranscriptFeature))
        transcript_attributes = re.split(r';',transcript_feature[8])

        transcript_attributes.insert(0,transcript_attributes[-3].lstrip(),) #for transcripts, gene_id is placed in the penultimate element [to access, antepenult, due to split ;]
        transcript_attributes[1] = " " + transcript_attributes[1] #add space in the transcript id field to be equal to all others
        del transcript_attributes[-3] #remove previous gene_id

        return transcript_feature,transcript_attributes


    def gene_id_rearrrangement_in_exons(GffutilsExonFeature):
        exon_feature = re.split(r'\t',str(children))
        exon_attributes = re.split(r';',exon_feature[8])

        exon_attributes.insert(0,exon_attributes[-2].lstrip(),) #for exons, gene_id is placed in the last element [to access, penultimante, due to split ;]
        exon_attributes[1] = " " + exon_attributes[1] #add space in the transcript id field to be equal to all others
        del exon_attributes[-2]

        return exon_feature,exon_attributes


    with open(outputbasename + 'merged.gtf', "w") as file:

        if software == "cufflinks" or software == "stringtie": #write transcript and exons
            for gene in db.features_of_type("gene"):

                for transcript in db.children(gene, level=1, featuretype='transcript' ):
                    feature,newAttribute = gene_id_rearrangement_in_transcripts(transcript)
                    file.write(newattribute(feature,newAttribute))

                    for children in db.children(transcript, level=1, featuretype='exon'):
                        feature,newAttribute = gene_id_rearrrangement_in_exons(children)
                        file.write(newattribute(feature,newAttribute))



        elif software == "cuffmerge":
            for transcript in db.features_of_type("transcript"):

                for children in db.children(transcript, level=1, featuretype='exon' ):
                    feature,newAttribute = gene_id_rearrrangement_in_exons(children)
                    file.write(newattribute(feature,newAttribute))


    file.close()

def generalStats(db):
    number_of_genes = db.count_features_of_type("gene")
    number_of_transcripts=db.count_features_of_type("transcript")
    number_of_exons=db.count_features_of_type("exon")

    return '{}'.format(number_of_genes),'{}'.format(number_of_transcripts),'{}'.format(number_of_exons)

def writeOutputStats(outputbasename,db,numb_repeated,original_numb_exon):
    if os.path.exists(outputbasename + 'stats.txt'):
        os.remove(outputbasename + 'stats.txt')

    a = generalStats(db)
    print type(a)
    #with open(outputbasename + '-info.txt', "w") as csvfile:
    #     writer = csv.writer(csvfile,dialect=csv.excel_tab)
    #     writer.writerow(('#Total umber of SNPs in differential expressed genes:', target_snps))
    #     writer.writerow(('#Number of genes with any SNP found', len(genes_with_snp)))
    #     writer.writerow('')
    #     writer.writerow(('#List of genes with SNPs:',''))
    #     for gene in genes_with_snp:
    #         writer.writerow((gene, ''))
    #     csvfile.close()

def transcriptExonUsage(db):

    fout = open('/home/pedro/Desktop/new.bed','w')
    exons= db.features_of_type("exon")

    for gene in db.features_of_type("gene"):

        print gene
        transcripts = db.children(gene.id, level=1, featuretype='transcript')
        for transcript in transcripts:
            #write bed to file and perform merge within python with subprocess module
#            fout.write(db.bed12(transcript)+"\n")
    #        print(transcript)
     #       for exon in db.children(transcript, featuretype='exon', order_by='start'):
      #          print(exon)
        #generate several stats according the tutorial
  #          merged_transcripts=db.merge(transcripts, ignore_strand=True)
   #         for exons in merged_transcripts:
             #  db.children(merged_exons, featuretype='exon', order_by='start'):
                print exons
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



def createGffUtilsCuffmerge(gtf_file,forceNewDB, isVerbose):
    dbname=os.path.basename(gtf_file).split('.')[0] + "DB"
    dialect=helpers.infer_dialect(['Potrx000002	Cufflinks	exon	8052	8625	.	-	.	gene_id "XLOC_000003"; transcript_id "TCONS_00000005"; exon_number "2"; '
                                   'gene_name "Potrx000002g00030"; oId "Potrx000002g00030.1"; nearest_ref "Potrx000002g00030.1"; class_code "="; tss_id "TSS3"; p_id "P3";'])


    try:

        db=gffutils.create_db(gtf_file, dbfn=dbname, id_spec={'gene': ['gene_id', 'gene_name', 'nearest_ref'],'transcript' : ['transcript_id', 'oId'], 'exon': 'exon_number'},
                              merge_strategy="merge",keep_order=True, sort_attribute_values=True, disable_infer_transcripts=False, disable_infer_genes=False,
                              dialect=dialect, checklines=500 ,verbose=isVerbose,force=forceNewDB)
    except:
        logging.warning("Database already exists. Will use the current one.")
        db=gffutils.FeatureDB(dbname, keep_order=True)

#    transcriptExonUsage(db)
    return db


def createGffUtilsCufflinks(gtf_file,forceNewDB, isVerbose):
    dbname=os.path.basename(gtf_file).split('.')[0] + "DB"
    dialect=helpers.infer_dialect(['Potrx000002	Cufflinks	exon	2615	2977	1	+	.	gene_id "Potrx000002g00010"; transcript_id "Potrx000002g00010.1"; exon_number "1"; '
                                   'FPKM "0.0000000000"; frac "0.000000"; conf_lo "0.000000"; conf_hi "0.000000"; cov "0.000000";'])

    try:
        db=gffutils.create_db(gtf_file, dbfn=dbname, id_spec={'gene': ['gene_id', 'gene_name', 'nearest_ref'],'transcript' : ['transcript_id', 'oId'], 'exon': 'exon_number'},
                              merge_strategy="merge",keep_order=True, sort_attribute_values=True ,disable_infer_transcripts=True, disable_infer_genes=False,
                              dialect=dialect, checklines=500 ,verbose=isVerbose,force=forceNewDB)
    except:
        logging.warning("Database already exists. Will use the current one.")
        db=gffutils.FeatureDB(dbname, keep_order=True)

    return db


def createGffUtilsStringtie(gtf_file,forceNewDB, isVerbose):
    dbname=os.path.basename(gtf_file).split('.')[0] + "DB"
    dialect=helpers.infer_dialect(['Potrx000002 StringTie	exon	8052	8625	1000	-	.	gene_id "STRG.2"; transcript_id "STRG.2.1";'
                                   ' exon_number "2"; reference_id "TCONS_00000004"; ref_gene_id "XLOC_000003"; ref_gene_name "Potrx000002g00030"; cov "2.764808";'])

    try:
        db=gffutils.create_db(gtf_file, dbfn=dbname, id_spec={'gene': ['gene_id', 'ref_gene_id', 'ref_gene_name'],'transcript' : ['transcript_id', 'reference_id'], 'exon': 'exon_number'},
                              merge_strategy="merge", keep_order=True, sort_attribute_values=True, disable_infer_transcripts=True, disable_infer_genes=False, dialect=dialect,
                              checklines=500 ,verbose=isVerbose,force=forceNewDB)
    except:
        logging.warning("Database already exists. Will use the current one.")
        db=gffutils.FeatureDB(dbname, keep_order=True)
    return db




def main():
    parser = argparse.ArgumentParser(description='Script to check the relationship between genes, transcripts and exons in a GTF file. Requires the gffutils module to be installed. ')
    parser.add_argument(dest='gtf_file', metavar='gtf_file', nargs=1, help='GTF file to be processed.')
    parser.add_argument(dest='software', metavar='software', nargs=1, type=str, choices=['cuffmerge','cufflinks','stringtie'],help='Tool that produced the input gtf file.')
    parser.add_argument(dest='output_prefix', metavar='output_prefix', nargs=1, help='Basename for the ouptut files.')
    parser.add_argument('--writeNewGTF', action='store_true', help='Write a new GTF file where the same exons are merged into the same feature.')
    parser.add_argument('--force', action='store_true', help='Overewrite existing gffutils database.')
    parser.add_argument('--verbose', action='store_true', help='Print more infomation while creating gffutils database.')

    args = parser.parse_args()



    if "cuffmerge" in args.software:

        db = createGffUtilsCuffmerge(args.gtf_file[0], args.force, args.verbose)
        numb_repeated, list_repeated,original_numb_exon = getRepeatedTranscripts(args.gtf_file[0])
        writeOutputStats(args.output_prefix[0],db,numb_repeated,original_numb_exon)
        if args.writeNewGTF:
            newGTFmerged(db,args.output_prefix[0],args.software[0])

    elif "cufflinks" in args.software:
        db = createGffUtilsCufflinks(args.gtf_file[0], args.force, args.verbose)
        numb_repeated, list_repeated,original_numb_exon = getRepeatedTranscripts(args.gtf_file[0])
        writeOutputStats(args.output_prefix[0],db,numb_repeated,original_numb_exon)
        if args.writeNewGTF:
            newGTFmerged(db,args.output_prefix[0],args.software[0])

    elif "stringtie" in args.software:
        db = createGffUtilsStringtie(args.gtf_file[0],args.force, args.verbose)
        numb_repeated, list_repeated,original_numb_exon = getRepeatedTranscripts(args.gtf_file[0])
        writeOutputStats(args.output_prefix[0],db,numb_repeated,original_numb_exon)
        if args.writeNewGTF:
            newGTFmerged(db,args.output_prefix[0],args.software[0])


if __name__ == "__main__":
    main()