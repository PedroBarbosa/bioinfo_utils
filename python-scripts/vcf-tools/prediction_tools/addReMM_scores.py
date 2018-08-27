import argparse
import os.path
from cyvcf2 import VCF, Writer
import pysam

def checkIndex(remmfile):
    if not os.path.isfile(remmfile + '.tbi'):
        print('Index file is required. Make sure that {} file exists in the same path'.format(remmfile + '.tbi'))
        exit(1)

def processVCF(invcf,remm,out):
    vcf_data = VCF(invcf, gts012=True)
    tbx = pysam.TabixFile(remm)
    vcf_data.add_info_to_header({'ID': 'ReMM', 'Description': 'The Regulatory Mendelian Mutation (ReMM) for relevance prediction of non-coding variations (SNVs and small InDels)','Type':'String', 'Number': '.'})
    w = Writer(out, vcf_data)
    for record in vcf_data:
        try:
            for row in tbx.fetch(record.CHROM.replace('chr',''), record.start, record.end):
                if int(str(row).split()[1]) == record.POS:
                    record.INFO["ReMM"] = str(row).split()[2]
            if not record.INFO["ReMM"]:
                record.INFO["ReMM"] = "."
        except ValueError:
            record.INFO["ReMM"] = "."
        w.write_record(record)
def main():
    parser = argparse.ArgumentParser(description='Script to append ReMM scores into a VCF')
    parser.add_argument(dest='vcf', help='Path to the vcf file')
    parser.add_argument(dest='remm', help='Path to the ReMM scores file. A tabix indexed file is required')
    parser.add_argument(dest='out',help='Path to the VCF output file')


    args = parser.parse_args()

    checkIndex(args.remm)
    processVCF(args.vcf,args.remm,args.out)

if __name__ == "__main__":
    main()
