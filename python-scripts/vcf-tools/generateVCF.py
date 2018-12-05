import argparse
from cyvcf2 import VCF, Writer
import numpy as np
import itertools
import random

def get_number_individuals(gnomad_vcf,pop):
    vcf_data = VCF(gnomad_vcf, gts012=True)
    nind = [int(record.INFO["AN"]/2) if pop == "All" else int(record.INFO["AN" + "_" + pop] / 2) for record in vcf_data]
    return int(np.median(nind))

def simulate_genotypes(nind,nhomalt,nhet):
    homalt_genotypes = [i for i in itertools.repeat([1,1,False], nhomalt)]
    het_genotypes = [i for i in itertools.repeat([0, 1, False], nhet)]
    ref_genotypes = [i for i in itertools.repeat([0, 0, False], nind - nhomalt - nhet)]
    genotypes = [l for l in [ref_genotypes,homalt_genotypes,het_genotypes] if l]
    l = list(itertools.chain(*genotypes))
    merged = random.sample(l, len(l))

def generate_vcf(gnomad_vcf,pop):
    nind = get_number_individuals(gnomad_vcf,pop)
    vcf_data = VCF(gnomad_vcf, gts012=True)
    for record in vcf_data:
        info=record.INFO
        if pop == "All":
            nhomalt = info["nhomalt"]
            nheterozygous = info["AC"] - (nhomalt * 2)
        else:
            nhomalt = info["nhomalt" + "_" + pop]
            nheterozygous = info["AC" + "_" + pop] - (nhomalt * 2)

        simulate_genotypes(nind,nhomalt,nheterozygous)


        #print(record.genotypes)#format("DP")))
        #print(nheterozygous)
        #print(nhomalt)



def main():
    parser = argparse.ArgumentParser(description='Script to create add genotypes into a VCF file.')
    parser.add_argument(dest='vcf_gnomad', help='Path to the gnomAD vcf file. This file must have been split and normalized before.')
    parser.add_argument(dest='vcf_absent_in_gnomad', help='Path to the vcf file containing the variants absent in gnmoAD after intersecting with the case-study vcf')
    parser.add_argument('-p','--population', default="All",choices=("All","raw","nfe","nfe_nwe","nfe_seu","amr","afr","eas"), help='Additional fields to include')
    args = parser.parse_args()

    generate_vcf(args.vcf_gnomad,args.population)

if __name__ == "__main__":
    main()