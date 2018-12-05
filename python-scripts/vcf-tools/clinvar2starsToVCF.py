import argparse
import pandas as pd
from cyvcf2 import VCF
from collections import Counter


def clinvarVCFReader(vcf,clinvarids,outvcf):
    with open(outvcf,'w') as out:
        [(out.write(str(variant)),clinvarids.remove(int(variant.ID))) for variant in VCF(vcf) if int(variant.ID) in clinvarids]
    out.close()
    with open("Invalid_dataset_IDs.txt",'w') as out:
        for item in clinvarids:
            out.write("{}\n".format(item))

def datasetReader(infile):
    df = pd.read_csv(infile,sep=",", low_memory=False)
    records=[]
    for index, row in df.iterrows():
        try:
            #records.append(("chr" + row["chrom"],int(row["hg19.start"]),row["clinvar.rsid"],row["clinvar.ref"],row["clinvar.alt"],".",".","."))
            records.append(row["clinvar.variant_id"])
        except ValueError:
            a=1
            #print(row["chrom"],row["clinvar.hg19.start"], row["clinvar.rsid"],row["clinvar.ref"], row["clinvar.alt"],".",".",".")
    return records

def main():
    parser = argparse.ArgumentParser(
        description='Script to generate VCF file from the alcides dataset')
    parser.add_argument(dest='dataset',  help='Path to the dataset file')
    parser.add_argument(dest='clinvarvcf', help='Path to the clinvar vcf file')
    parser.add_argument(dest='outvcf', help='Path to the output vcf file')
    args = parser.parse_args()

    clinvarids=datasetReader(args.dataset)
    clinvarVCFReader(args.clinvarvcf,clinvarids,args.outvcf)
if __name__ == "__main__":
    main()