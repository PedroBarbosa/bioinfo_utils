import argparse
from cyvcf2 import VCF
from collections import defaultdict
def parseVCF(invcf,outbasename,reportMultipleSamples):
    vcf_data = VCF(invcf,gts012=True)
    samples=vcf_data.samples
    multiple_samples=defaultdict(list)
    with open(outbasename + "_genotypes.txt",'w') as out:
        out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("CHROM","POS","REF","ALT","NUM_HET","HET_SAMPLES","NUM_HOM_ALT","HOM_ALT_SAMPLES"))
        for record in vcf_data:
            home_var= (record.gt_types==2).nonzero()[0]
            heterozygous = (record.gt_types==1).nonzero()[0]
            samples_het,samples_homvar = [],[]
            if heterozygous.size:
                [(samples_het.append(samples[i]),multiple_samples[samples[i]].append((record.CHROM,str(record.POS)))) for i in heterozygous]
            if home_var.size:
                [(samples_homvar.append(samples[i]),multiple_samples[samples[i]].append((record.CHROM,str(record.POS)))) for i in home_var]
            out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(record.CHROM, record.POS, record.REF, record.ALT[0],
                                                        record.num_het, ';'.join(samples_het),record.num_hom_alt,';'.join(samples_homvar)))
    out.close()
    if multiple_samples and reportMultipleSamples:
        outlist = [[k,v]for k, v in multiple_samples.items() if len(v) > 1]
        if not outlist:
            print("{}".format("No samples carry more than one variant represented in the input VCF."))
        else:
            with open(outbasename + "_multipleVariantSamples.txt",'w') as outmult:
                outmult.write("{}\t{}\t{}\n".format("Sample","#variants","variant_id"))
                for sample in outlist:
                    outmult.write("{}\t{}\t{}\n".format(sample[0],len(sample[1]),';'.join('_'.join(i) for i in sample[1])))
            outmult.close()



def main():
    parser = argparse.ArgumentParser(description='Script to process genotype VCF fields in a sample-based manner')
    parser.add_argument(dest='vcf', help='Path to the vcf')
    parser.add_argument(dest='outbasename', help='Output basename')
    parser.add_argument('-m','--multipleSamples', action='store_true', help='Flag to report samples that carry multiple variants')
    args = parser.parse_args()

    parseVCF(args.vcf,args.outbasename,args.multipleSamples)
if __name__ == "__main__":
    main()