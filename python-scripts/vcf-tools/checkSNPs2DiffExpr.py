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




def processVCF(vcf,geneList):
    target_snps=0
    genes_with_snp=set()
    with open(vcf) as file:

        for line in file:

            if not line.startswith("#"):
                contig_name = line.split()[0]
                if contig_name in geneList:

                    genes_with_snp.add(contig_name)
                    target_snps += 1


    return genes_with_snp,target_snps


def writeOutput(genes_with_snp,target_snps,outputfile):

    if os.path.exists(outputfile):
        os.remove(outputfile)
    with open(outputfile, "w") as csvfile:
        writer = csv.writer(csvfile,dialect=csv.excel_tab)
        writer.writerow(('#Total umber of SNPs in differential expressed genes:', target_snps))
        writer.writerow(('#Number of genes with any SNP found', len(genes_with_snp)))
        writer.writerow('')
        writer.writerow(('#List of genes with SNPs:',''))
        for gene in genes_with_snp:
            writer.writerow((gene, ''))

    csvfile.close()


parser = argparse.ArgumentParser(description='Script to check over a list of differential expressed genes the existence of SNPs.')
parser.add_argument(dest='vcf_file', metavar='vcf_file', nargs=1,help='File to be processed.')
parser.add_argument(dest='diffExpr_genes', metavar='genes_file', nargs=1, help='File with the list of genes up and down regulated')
parser.add_argument(dest='output_file', metavar='output_file', nargs=1, help='Output file')
parser.add_argument('-t', '--isTranscoder', action='store_true', help='Parse gene names if gene predictions were made with transcoder software')
args = parser.parse_args()

if __name__ == "__main__":

    if args.isTranscoder:
        geneList=getGeneListFromTranscoder(args.diffExpr_genes[0])
    else:
        geneList=getGeneList(args.diffExpr_genes[0])
    print("Scanning vcf file..")
    genes_with_snp,target_snps = processVCF(args.vcf_file[0], geneList)
    print("Writing output..")
    writeOutput(genes_with_snp,target_snps,args.output_file[0])
    print("Done.")
