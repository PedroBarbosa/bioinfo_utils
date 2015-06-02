__author__ = 'pedro'

import argparse
import os
import csv


def getGeneListFromTranscoder(file):
    geneList=[]
    with open(file) as file:
        for line in file:
	   if not line.startswith("#"):
           	geneList.append(line.split("|")[0].rstrip())

    return geneList




def getGeneList(file):
    geneList=[]
    with open(file) as file:
        for line in file:
            if not line.startswith("#"):
		geneList.append(line.rstrip())
    return geneList




def processAnnotatedVCF(vcf,geneList):

    final_data = []
    with open(vcf) as file:

        for line in file:

            if not line.startswith("#"):
                transcript_name = line.split()[0]
                transcript_snp_position = line.split()[1]
                info_fields = line.split()[7]
                if "ANN=" in info_fields and transcript_name in geneList:
                    CDSs_annotation = info_fields.split(",")
                    for CDS in CDSs_annotation:
                        effect = CDS.split("|")[1]
                        impact = CDS.split("|")[2]
                        final_data.append(transcript_name + ";yes;" + transcript_snp_position + ";" + effect + ";" + impact)

    return final_data


def writeOutput(final_data, geneList,outputprefix):

    genes_processed = set()
    ##info file
    if os.path.exists(outputprefix+'.tsv'):
        os.remove(outputprefix+'.tsv')
    with open(outputprefix+'.tsv', "w") as csvfile:
        writer = csv.writer(csvfile,dialect=csv.excel_tab)
        writer.writerow(('#Transcript_id', 'has_SNP?', 'transcript_position', 'SNP_effect', 'SNP_impact'))
        for snp in final_data:
            info = snp.split(';')
            writer.writerow((info[0],info[1],info[2],info[3],info[4]))
            genes_processed.add(info[0])

        #add the genes with no snps found
        for gene in geneList:
            if gene not in genes_processed:
                writer.writerow((gene, "no", "-", "-", "-"))
                
        csvfile.close()






parser = argparse.ArgumentParser(description='Script to create a table with information regarding SNPs for all the differential expressed genes in a study.')
parser.add_argument(dest='diffExpr_genes', metavar='genes_file', nargs=1, help='File with the list of genes up and down regulated.')
parser.add_argument(dest='annotated_vcf_file', metavar='ann_vcf_file', nargs=1,help='Variant annotated file.')
parser.add_argument(dest='output_prefix', metavar='output_prefix', nargs=1, help='Basename to write on the ouptut files.')
parser.add_argument('-t', '--isTransDecoder', action='store_true', help='Parse gene names if gene predictions were made with TransdDecoder software.')
args = parser.parse_args()

if __name__ == "__main__":

    if args.isTransDecoder:
        geneList=getGeneListFromTranscoder(args.diffExpr_genes[0])
    else:
        geneList=getGeneList(args.diffExpr_genes[0])
    print("Scanning vcf file..")
    final_data = processAnnotatedVCF(args.annotated_vcf_file[0], geneList)
    print("Writing output..")
    writeOutput(final_data, geneList, args.output_prefix[0])
    print("Done.")

