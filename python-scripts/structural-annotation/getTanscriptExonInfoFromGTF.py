__author__ = 'pedro'

import argparse
import os
import gffutils
from gffutils import helpers
import logging
import re
import csv
import sys
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')

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
    logging.info("Checking for repeated transcripts..[Deprecated function]")
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
                    logging.error("ERROR with input file.Please check your gtf file, might exist some empty lines at the end of the file.\n")
                    sys.exit(1)

        if exons_string in exonsPertranscript: #Last line of file
            numb_repeated_transcripts += 1
            repeated_transcripts.append(transcript_id) #not previous because file is read


    return numb_repeated_transcripts, repeated_transcripts,original_numb_exon



def newGTFmerged(db,outputbasename,software):
    logging.info("Creating new GTF file with the merged exons..")
    if os.path.exists(outputbasename + '-merged.gtf'):
        os.remove(outputbasename + '-merged.gtf')

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


    with open(outputbasename + '-merged.gtf', "w") as file:

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

    logging.info("Calculating general stats:")
    def exonsAverage(db):
        logging.info("\tPer feature exon average..")
        gene_exon_count = 0
        gene_count = 0
        transcript_count = 0
        transcript_exon_count = 0
        for gene in db.features_of_type('gene'):
            exonsInGene = 0

            for transcript in db.children(gene,1):
                exonsInTranscript = 0
                if transcript.featuretype == 'transcript':

                    for child in db.children(transcript,1):
                        if child.featuretype == 'exon':
                            exonsInTranscript += 1


            transcript_exon_count += exonsInTranscript
            transcript_count += 1


            # get all grandchildren, only counting the exons
            for child in db.children(gene.id,2):
                if child.featuretype == 'exon':
                    exonsInGene += 1

            gene_exon_count += exonsInGene
            gene_count += 1

        meanExonGene = round(float(gene_exon_count) / gene_count,2)
        meanExonTranscript = round(float(transcript_exon_count)/ transcript_count,2)

        return meanExonGene, meanExonTranscript

    def constitutive_exons(db):
        logging.info("\tConstitutive exons..")
        constitutive_exons = []
        exone=0
        for exon in db.features_of_type('exon'):

             geneID = re.findall(r'\'([^\']*)\'', str(exon['gene_id']))

             n_iso =len(list(db.children(geneID[0], featuretype= 'transcript')))
             exone += 1

             if len(list(db.parents(exon, featuretype='transcript'))) == n_iso:
                constitutive_exons.append(str(exon.id))

        # for gene in db.features_of_type('gene'):
        #     exone +=1
        #     n_isoforms = len(list(db.children(gene, featuretype='transcript')))
        #
        #     for exon in db.children(gene, 2, featuretype='exon'):
        #         parents = db.parents(exon, featuretype='transcript')
        #
        #         if len(list(parents)) == n_isoforms:
        #             constitutive_exons.append(str(exon.id))

        return constitutive_exons

    number_of_genes = str(db.count_features_of_type("gene"))
    number_of_exons=str(db.count_features_of_type("exon"))
    number_of_transcripts= str(db.count_features_of_type("transcript"))



    logging.info("\tAverage features length..")
    gene_lengths, max_gene_len, total_transcript_length, total_exons_length = 0,0,0,0

    for gene in db.features_of_type('gene'):
        if len(gene) > max_gene_len:
            max_gene_len = len(gene)
            longest_gene = gene.id
        gene_lengths += len(gene)

        transcript_lengths = 0
        exons_lengths = 0
        for transcript in db.children(gene,featuretype='transcript'):
            transcript_lengths += len(transcript)

        ## add per transcript average per gene and then divide per the total number of genes. Intuitive average would show greater transcript average than gene average length.
        total_transcript_length += round(float(transcript_lengths) / int(len(list(db.children(gene, featuretype='transcript')))),2)


        for exon in db.children(gene, featuretype='exon'):
            exons_lengths += len(exon)
        total_exons_length += round(float(exons_lengths) / int(len(list(db.children(gene, featuretype='exon')))),2)


    mean_gene_length = round(float(gene_lengths) / int(number_of_genes),2)
    mean_transcript_length = round(float(total_transcript_length) / int(number_of_genes),2)
    mean_exon_length = round(float(total_exons_length) / int(number_of_genes),2)


    avg_transcriptPerGene = round(float(number_of_transcripts)/int(number_of_genes),2)
    meanExonGene, meanExonTranscript = exonsAverage(db)
    constitutive_exons = constitutive_exons(db)

    return number_of_genes, number_of_transcripts, number_of_exons, mean_gene_length, mean_transcript_length, mean_exon_length, max_gene_len, longest_gene, meanExonGene,avg_transcriptPerGene, meanExonTranscript, constitutive_exons
    #return number_of_genes, number_of_transcripts, number_of_exons, mean_gene_length, mean_transcript_length, mean_exon_length, max_gene_len, longest_gene, meanExonGene,avg_transcriptPerGene, meanExonTranscript
    #return number_of_genes, number_of_transcripts, number_of_exons, mean_gene_length, mean_transcript_length, mean_exon_length, max_gene_len, longest_gene,avg_transcriptPerGene


def writeOutputStats(outputbasename,db,numb_repeated,original_numb_exon):


    if os.path.exists(outputbasename + '-info.txt'):
        os.remove(outputbasename + '-info.txt')

    mainStats = generalStats(db)

    logging.info("Writing stats to file..")
    with open(outputbasename + '-info.txt', "w") as csvfile:
         writer = csv.writer(csvfile,dialect=csv.excel_tab)
         writer.writerow(('#Number of genes in GTF:', mainStats[0]))
         writer.writerow(('#Number of transcripts in GTF:', mainStats[1]))
         writer.writerow(('#Number of exons in original GTF:', original_numb_exon))
         writer.writerow(('#Number of exons after merging the one with the same coordinates:', mainStats[2]))
         writer.writerow('')
         writer.writerow(('#Longest gene [length]', mainStats[7] + '[' + str(mainStats[6]) + ']'))
         writer.writerow(('#Average gene length', mainStats[3]))
         writer.writerow(('#Average transcript length [with introns spanning corresponding exons included]', mainStats[4]))
         writer.writerow(('#Average exon length', mainStats[5]))
         writer.writerow('')
         writer.writerow(('#Average number of exons per gene', mainStats[8]))
         writer.writerow(('#Average number of transcripts per gene', mainStats[9]))
         writer.writerow(('#Average number of exons per expressed transcript', mainStats[10]))
         writer.writerow(('#Number of constitutive exons [present in all isoforms of a gene]', round(len(mainStats[11]),4)))
         writer.writerow('')
         writer.writerow(('#Number of repeated transcripts in GTF with the same string of exons [does not mean they are the same][Deprecated]:',numb_repeated))

         csvfile.close()



def createGffUtilsCuffmerge(gtf_file,forceNewDB, isVerbose):
    dbname=os.path.abspath(gtf_file).split('.')[0] + "DB"
    dialect=helpers.infer_dialect(['Potrx000002	Cufflinks	exon	8052	8625	.	-	.	gene_id "XLOC_000003"; transcript_id "TCONS_00000005"; exon_number "2"; '
                                   'gene_name "Potrx000002g00030"; oId "Potrx000002g00030.1"; nearest_ref "Potrx000002g00030.1"; class_code "="; tss_id "TSS3"; p_id "P3";'])

    try:

        db=gffutils.create_db(gtf_file, dbfn=dbname, id_spec={'gene': ['gene_id', 'gene_name', 'nearest_ref'],'transcript' : ['transcript_id', 'oId'], 'exon': 'exon_number'},
                              merge_strategy="merge",keep_order=True, sort_attribute_values=True, disable_infer_transcripts=False, disable_infer_genes=False,
                              dialect=dialect, checklines=500 ,verbose=isVerbose,force=forceNewDB)
        logging.info("Database " + dbname + "successfuly generated.")
    except:
        logging.info("Database already exists. Will use the current one.")
        try:
            db=gffutils.FeatureDB(dbname, keep_order=True)
        except TypeError:
            logging.error("Previous generated database might be corrupted. Please set --force to overwrite the database and --verbose if you want to follow the steps of gffutils"
                          " database creation.")
            exit(1)
    return db


def createGffUtilsCufflinks(gtf_file,forceNewDB, isVerbose):
    dbname=os.path.abspath(gtf_file).split('.')[0] + "DB"
    dialect=helpers.infer_dialect(['Potrx000002	Cufflinks	exon	2615	2977	1	+	.	gene_id "Potrx000002g00010"; transcript_id "Potrx000002g00010.1"; exon_number "1"; '
                                   'FPKM "0.0000000000"; frac "0.000000"; conf_lo "0.000000"; conf_hi "0.000000"; cov "0.000000";'])

    try:
        db=gffutils.create_db(gtf_file, dbfn=dbname, id_spec={'gene': ['gene_id'],'transcript' : ['transcript_id'], 'exon': 'exon_number'},
                              merge_strategy="merge",keep_order=True, sort_attribute_values=True ,disable_infer_transcripts=True, disable_infer_genes=False,
                              dialect=dialect, checklines=500 ,verbose=isVerbose,force=forceNewDB)

        logging.info("Database " + dbname + "successfuly generated.")
    except:
        logging.info("Database already exists. Will use the current one.")
        try:
            db=gffutils.FeatureDB(dbname, keep_order=True)
        except TypeError:
            logging.error("Previous generated database might be corrupted. Please set --force to overwrite the database and --verbose if you want to follow the steps of gffutils"
                          " database creation.")
            exit(1)


    return db


def createGffUtilsStringtie(gtf_file,forceNewDB, isVerbose):
    dbname=os.path.abspath(gtf_file).split('.')[0] + "DB"
    dialect=helpers.infer_dialect(['Potrx000002 StringTie	exon	8052	8625	1000	-	.	gene_id "STRG.2"; transcript_id "STRG.2.1";'
                                   ' exon_number "2"; reference_id "TCONS_00000004"; ref_gene_id "XLOC_000003"; ref_gene_name "Potrx000002g00030"; cov "2.764808";'])

    try:
        db=gffutils.create_db(gtf_file, dbfn=dbname, id_spec={'gene': ['gene_id', 'ref_gene_id', 'ref_gene_name'],'transcript' : ['transcript_id', 'reference_id'], 'exon': 'exon_number'},
                              merge_strategy="merge", keep_order=True, sort_attribute_values=True, disable_infer_transcripts=True, disable_infer_genes=False, dialect=dialect,
                              checklines=500 ,verbose=isVerbose,force=forceNewDB)
        logging.info("Database " + dbname + "successfuly generated.")

    except:
        logging.info("Database already exists. Will use the current one.")
        try:
            db=gffutils.FeatureDB(dbname, keep_order=True)
        except TypeError:
            logging.error("Previous generated database might be corrupted. Please set --force to overwrite the database and --verbose if you want to follow the steps of gffutils"
                          " database creation.")
            exit(1)

    return db




def main():
    parser = argparse.ArgumentParser(description='Script to check and generate stats about the relationship between genes, transcripts and exons in a GTF file. Requires the gffutils module to be installed. ')
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