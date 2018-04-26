import argparse
from cyvcf2 import VCF

def parseVCF(invcf,outfile):
    vcf_data = VCF(invcf,gts012=True)
    samples=vcf_data.samples
    with open(outfile,'w') as out:
        out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("CHROM","POS","REF","ALT","NUM_HET","HET_SAMPLES","NUM_HOM_ALT","HOM_ALT_SAMPLES"))
        for record in vcf_data:
            home_var= (record.gt_types==2).nonzero()[0]
            heterozygous = (record.gt_types==1).nonzero()[0]
            samples_het,samples_homvar = [],[]
            if heterozygous.size:
                [samples_het.append(samples[i]) for i in heterozygous]
            if home_var.size:
                [samples_homvar.append(samples[i]) for i in home_var]
            out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(record.CHROM, record.POS, record.REF, record.ALT[0],
                                                        record.num_het, ';'.join(samples_het),record.num_hom_alt,';'.join(samples_homvar)))
def main():
    parser = argparse.ArgumentParser(description='Script to process genotype VCF fields in a sample-based manner')
    parser.add_argument(dest='vcf', help='Path to the vcf')
    parser.add_argument(dest='outfile', help='Path to the output file')
    args = parser.parse_args()

    parseVCF(args.vcf,args.outfile)
if __name__ == "__main__":
    main()