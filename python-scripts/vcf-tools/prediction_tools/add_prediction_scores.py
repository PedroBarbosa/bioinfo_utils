import argparse
import os.path
from cyvcf2 import VCF, Writer
import pysam

def checkIndex(remmfile):
    if not os.path.isfile(remmfile + '.tbi'):
        print('Index file is required. Make sure that {} file exists in the same path'.format(remmfile + '.tbi'))
        exit(1)

def processVCF(invcf,remm,dann,out):
    vcf_data = VCF(invcf, gts012=True)
    tbx_remm = pysam.TabixFile(remm)
    tbx_dann = pysam.TabixFile(dann)
    vcf_data.add_info_to_header({'ID': 'ReMM', 'Description': 'The Regulatory Mendelian Mutation (ReMM) for relevance prediction of non-coding variations (SNVs and small InDels)','Type':'String', 'Number': '.'})
    vcf_data.add_info_to_header({'ID': 'DANN', 'Description': 'A deep neural network aimed to recognize pathogenic variants by annotating genetic variants, especially in noncoding regions.','Type': 'String', 'Number': '.'})
    w = Writer(out, vcf_data)
    for record in vcf_data:
        try:
            for row in tbx_remm.fetch(record.CHROM, record.start, record.end):

                if int(str(row).split()[1]) == record.POS:
                    record.INFO["ReMM"] = str(row).split()[2]
            if not record.INFO["ReMM"]:
                record.INFO["ReMM"] = "."
        except ValueError:
            record.INFO["ReMM"] = "."

        try:
            for row in tbx_dann.fetch(record.CHROM, record.start, record.end):
                if int(row.split()[1]) == record.POS and row.split()[2] == record.REF and row.split()[3] == record.ALT[0]:     
                    record.INFO["DANN"] = round(float(row.split()[4]),3)
                    break   
                else:
                    record.INFO["DANN"] = "."
        except ValueError:
            record.INFO["DANN"] = "."

        w.write_record(record)
def main():
    parser = argparse.ArgumentParser(description='Script to append ReMM scores into a VCF')
    parser.add_argument(dest='vcf', help='Path to the vcf file')
    parser.add_argument(dest='out',help='Path to the VCF output file')
    parser.add_argument('--remm', default="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/reMM_genomiser/ReMM.v0.3.1.txt.bgz",help='Path to the ReMM scores file. A tabix indexed file is required')
    parser.add_argument('--dann', default="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/wgsa_dbNSFP/resources/DANN/DANN_hg19.txt.gz",help='Path to the DANN scores file. A tabix indexed file is required')

    args = parser.parse_args()
    checkIndex(args.remm)
    checkIndex(args.dann)
    processVCF(args.vcf,args.remm,args.dann,args.out)

if __name__ == "__main__":
    main()

