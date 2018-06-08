import argparse
from cyvcf2 import VCF
from collections import defaultdict

def extractFields(vcf,fields):
    indexes=defaultdict(list)
    vcf_data = VCF(vcf, gts012=True)
    samples = vcf_data.samples
    no_ref_genotypes=(1,2)
    sample_based=defaultdict(list)
    for record in vcf_data:
        heterozygous = (record.gt_types == 1).nonzero()[0].tolist()
        hom_alt = (record.gt_types == 2).nonzero()[0].tolist()
        merged=heterozygous+hom_alt
        for i in merged:
            sample_based[samples[i]].extend(record)



def main():
    parser = argparse.ArgumentParser(description='Script to produces excel readable file for annotated variants')
    parser.add_argument(dest='vcf', help='Path to the vcf')
    parser.add_argument('-f','--fields', nargs='+', help='Additional fields to include')
    args = parser.parse_args()

    extractFields(args.vcf,args.fields)
if __name__ == "__main__":
    main()
