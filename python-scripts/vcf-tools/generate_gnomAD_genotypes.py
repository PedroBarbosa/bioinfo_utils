import argparse
from cyvcf2 import VCF
import numpy as np
import itertools
import random
import copy
import Bio


def get_format_fields_from_vcf(vcf):
    vcf_data = VCF(vcf, gts012=True)
    for record in vcf_data:
        return record.FORMAT


def generate_putative_GQ_DP(fields, nind):
    if "GT" and "GQ" in fields:
        gt_quals = np.full(shape=nind, fill_value=50, dtype=np.int)
        gt_dp = np.full(shape=nind, fill_value=100, dtype=np.int)
        return gt_dp, gt_quals
    else:
        print("GQ and DP fields are not present in your non-gnomAD VCF file ? Please check the FORMAT field.")
        exit(1)


def get_number_individuals(gnomad_vcf, pop):
    vcf_data = VCF(gnomad_vcf, gts012=True)
    nind = [int(record.INFO["AN"] / 2) if pop == "All" else int(record.INFO["AN" + "_" + pop] / 2) for record in
            vcf_data]
    return int(np.max(nind))


def simulate_genotypes(nind, nhomalt, nhet):
    homalt_genotypes = [i for i in itertools.repeat("1/1", nhomalt)]
    het_genotypes = [i for i in itertools.repeat("0/1", nhet)]
    ref_genotypes = [i for i in itertools.repeat("0/0", nind - nhomalt - nhet)]
    genotypes = [l for l in [ref_genotypes, homalt_genotypes, het_genotypes] if l]
    l = list(itertools.chain(*genotypes))
    return random.sample(l, len(l))


def generate_vcf(gnomad_vcf, outfile, pop, format_fields):
    nind = get_number_individuals(gnomad_vcf, pop)
    gt_dp, gt_qual = generate_putative_GQ_DP(format_fields, nind)
    vcf_data = VCF(gnomad_vcf, gts012=True)
    if outfile.endswith(".gz"):
        output = Bio.bgzf.BgzfWriter(filename=outfile)
    else:
        output=outfile
    out = open(output,'w')

    vcf_data.add_format_to_header({'ID': 'GT', 'Description': 'Genotype', 'Type': 'String', 'Number': 1})
    vcf_data.add_format_to_header(
        {'ID': 'AD', 'Description': 'Allelic depths for the ref and alt alleles in the order listed',
         'Type': 'Integer', 'Number': 1})
    vcf_data.add_format_to_header(
        {'ID': 'DP', 'Description': 'Approximate read depth', 'Type': 'Integer', 'Number': 1})
    vcf_data.add_format_to_header({'ID': 'GQ', 'Description': 'Genotyp Quality', 'Type': 'Integer', 'Number': 1})
    vcf_data.add_format_to_header({'ID': 'PL',
                                   'Description': 'normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification',
                                   'Type': 'Integer', 'Number': 1})

    individuals = ["ind_" + str(i) for i in range(1, nind, 1)]
    header = filter(None, vcf_data.raw_header.split("\n"))
    final_header = [line + '\t' + '\t'.join(individuals) if line.startswith("#CHROM") else line for line in header]
    out.write('\n'.join(final_header) + "\n")
    info_fields = [field["ID"] for field in vcf_data.header_iter() if field["HeaderType"] == "INFO"]

    for record in vcf_data:
        print("here")
        info = record.INFO
        if pop == "All":
            nhomalt = info["nhomalt"]
            nheterozygous = info["AC"] - (nhomalt * 2)
        else:
            nhomalt = info["nhomalt" + "_" + pop]
            nheterozygous = info["AC" + "_" + pop] - (nhomalt * 2)

        record_combined_gt = simulate_genotypes(nind, nhomalt, nheterozygous)
        gt_phred_ll_homref = np.zeros((nind,), dtype=np.int)
        gt_phred_ll_het = np.full(shape=nind, fill_value=1500, dtype=np.int)
        gt_phred_ll_homalt = np.full(shape=nind, fill_value=1500, dtype=np.int)
        gt_alt_depth = np.zeros((nind,), dtype=np.int)
        gt_ref_depth = copy.deepcopy(gt_dp)
        fmt = []
        for i, gt in enumerate(record_combined_gt):
            if gt == "0/1":  # if het
                gt_phred_ll_het[i] = 0
                gt_phred_ll_homref[i] = 1500
                gt_alt_depth[i] = gt_ref_depth[i] = 50
            elif gt == "1/1":  # if homalt
                gt_phred_ll_homalt[i] = 1500
                gt_phred_ll_homref[i] = 1500
                gt_alt_depth[i] = 100
                gt_ref_depth[i] = 0

            fmt.append("{}:{},{}:{}:{}:{},{},{}".format(gt, gt_ref_depth[i], gt_alt_depth[i], gt_dp[i], gt_qual[i],
                                                        gt_phred_ll_homref[i], gt_phred_ll_het[i],
                                                        gt_phred_ll_homalt[i]))
        str_info = []
        for i in info_fields:
            try:
                str_info.append(i + "=" + str(record.INFO[i]))
            except KeyError as abs_fied:
                #print("Field {} absent".format(abs_fied))
                continue

        write_record = ['.' if v is None else v for v in
                        [record.CHROM, str(record.POS), record.ID, record.REF, record.ALT[0], str(record.QUAL),
                         record.FILTER,
                         ';'.join(str_info), "GT:AD:DP:GQ:PL"]]

        out.write('\t'.join(write_record + fmt) + "\n")

    vcf_data.close()
    out.close()
    return nind

def add_absent_records(vcf_absent_gnomad,outfile,nind):
    #format_fields = get_format_fields_from_vcf(vcf_absent_gnomad)
    gt, gt_dp, gt_ref_depth, gt_alt_depth,gt_qual = "0/0",100,100,0,50
    gt_phred_ll_homref,gt_phred_ll_het,gt_phred_ll_homalt = 0,1500,1500
    fmt = ["{}:{},{}:{}:{}:{},{},{}".format(gt, gt_ref_depth, gt_alt_depth, gt_dp, gt_qual,
                                            gt_phred_ll_homref, gt_phred_ll_het, gt_phred_ll_homalt)] * nind
    vcf_data = VCF(vcf_absent_gnomad, gts012=True)
    info_fields = [field["ID"] for field in vcf_data.header_iter() if field["HeaderType"] == "INFO"]
    with open(outfile,'a') as out:
        for record in vcf_data:
            str_info = []
            for i in info_fields:
                try:
                    str_info.append(i + "=" + str(record.INFO[i]))
                except KeyError:
                    continue
            write_record = ['.' if v is None else v for v in
                            [record.CHROM, str(record.POS), record.ID, record.REF, record.ALT[0], str(record.QUAL),
                             record.FILTER,
                             ';'.join(str_info), "GT:AD:DP:GQ:PL"]]

            out.write('\t'.join(write_record + fmt) + "\n")
        vcf_data.close()
    out.close()

def main():
    parser = argparse.ArgumentParser(description='Script to create add genotypes into a VCF file.')
    parser.add_argument(dest='vcf_gnomad',
                        help='Path to the gnomAD vcf file. This file must have been split and normalized before.')
    parser.add_argument(dest='vcf_absent_in_gnomad',
                        help='Path to the vcf file containing the variants absent in gnmoAD after intersecting with the case-study vcf')
    parser.add_argument(dest='output_file', help='Output file name')
    parser.add_argument('-p', '--population', default="All",
                        choices=("All", "raw", "nfe", "nfe_nwe", "nfe_seu", "amr", "afr", "eas"),
                        help='Additional fields to include')
    args = parser.parse_args()

    fields = get_format_fields_from_vcf(args.vcf_absent_in_gnomad)
    nindividuals = generate_vcf(args.vcf_gnomad, args.output_file, args.population, fields)
    add_absent_records(args.vcf_absent_in_gnomad, args.output_file, nindividuals)

if __name__ == "__main__":
    main()
