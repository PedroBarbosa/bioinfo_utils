import argparse
from cyvcf2 import VCF
from collections import defaultdict
import hgvs.parser
import hgvs.enums
import hgvs.exceptions


def get_location(x,hp):
    try:
        v = hp.parse_hgvs_variant(x.split(" ")[0])
        if v.type == "m":
            return 'mithocondrial'
        elif v.type == "n":
            return 'noncodingRNA'
        elif v.posedit.pos.start.base < 0:
            return '5primeUTR'
        elif v.posedit.pos.start.datum == hgvs.enums.Datum.CDS_END and v.posedit.pos.start.base > 0:
            return '3primeUTR'
        elif v.posedit.pos.start.base > 0 and v.posedit.pos.start.offset == 0:
            return 'coding' #coding
        elif v.posedit.pos.start.base > 0 and abs(v.posedit.pos.start.offset) >= 100:
            return 'deepintronic (>100bp)'
        else:
            return 'intronic (10to100bp)'

    except hgvs.exceptions.HGVSParseError:
        return 'likely_promotor_or_intergenic'


def get_annotations(record, tools, samples, sample_based, ind, hp, specific_fields=None):
    if specific_fields is None:
        var = [record.CHROM, str(record.POS), record.REF, record.ALT[0],
           record.INFO["ANN"].split(",")[0].split("|")[tools.index('Existing_variation')],
           record.INFO["ANN"].split(",")[0].split("|")[tools.index('MAX_AF')],
           record.INFO["ANN"].split(",")[0].split("|")[tools.index('MAX_AF_POPS')],
           record.INFO["ANN"].split(",")[0].split("|")[tools.index('AF')],
           record.INFO["ANN"].split(",")[0].split("|")[tools.index('gnomADg_AF')],
           record.INFO["ANN"].split(",")[0].split("|")[tools.index('gnomADg_AF_nfe')],
           record.INFO["ANN"].split(",")[0].split("|")[tools.index('gnomADg_AF_nfe_nwe')],
           # record.INFO["ANN"].split(",")[0].split("|")[tools.index('gnomADg_AF_afr')],
           # record.INFO["ANN"].split(",")[0].split("|")[tools.index('gnomADg_AF_eas')],
           record.INFO["GERP"],
           record.INFO["phastCons"],
           record.INFO["phyloP"],
           record.INFO["ANN"].split(",")[0].split("|")[tools.index('CADD_PHRED')],
           max(record.INFO["ANN"].split(",")[0].split("|")[tools.index('ada_score')],
               record.INFO["ANN"].split(",")[0].split("|")[tools.index('rf_score')]),
           record.INFO["ANN"].split(",")[0].split("|")[tools.index('MaxEntScan_diff')]]

        for cons in record.INFO["ANN"].split(","):
            fields = cons.split("|")
            if len(fields) > 1:
                var.extend((fields[tools.index('Consequence')], fields[tools.index('IMPACT')], fields[tools.index('SYMBOL')],
                        fields[tools.index('Feature')],
                        fields[tools.index('HGVSg')],
                        fields[tools.index('HGVSc')],
                        fields[tools.index('HGVSp')],
                        get_location(fields[tools.index('HGVSc')], hp),
                        fields[tools.index('CLIN_SIG')],
                        ))

        sample_based[samples[ind]].append(var)
        return sample_based

    else:
        var = [record.CHROM, str(record.POS), record.REF, record.ALT[0],
               record.INFO["ANN"].split(",")[0].split("|")[tools.index('Existing_variation')],
               record.INFO["ANN"].split(",")[0].split("|")[tools.index('MAX_AF')],
               record.INFO["ANN"].split(",")[0].split("|")[tools.index('MAX_AF_POPS')],
               record.INFO["ANN"].split(",")[0].split("|")[tools.index('AF')],
               record.INFO["ANN"].split(",")[0].split("|")[tools.index('gnomADg_AF')],
               record.INFO["ANN"].split(",")[0].split("|")[tools.index('gnomADg_AF_nfe')],
               record.INFO["ANN"].split(",")[0].split("|")[tools.index('gnomADg_AF_nfe_nwe')],
               get_location(record.INFO["ANN"].split(",")[0].split("|")[tools.index('HGVSc')], hp)]

        for field in specific_fields:
            res=""
            try:
                res = record.INFO["ANN"].split(",")[0].split("|")[tools.index(field)]
            except ValueError:
                try:
                    res = record.INFO[field]
                except KeyError:
                    print("Record {},{},{},{} does not have {} in the INFO or ANN field/subfield, respectively."
                          .format(record.CHROM, str(record.POS), record.REF, record.ALT[0], field))
                    exit(1)

            var.extend([res])

        sample_based[samples[ind]].append(var)
        return sample_based


def write_output(sample_based_dict, output_file, negative_samples_provided, specific_fields=None):
    with open(output_file, 'w') as out:
        if specific_fields is None:
            s = ["Sample_ID", "isPositive"] if negative_samples_provided else ["Sample_ID"]
            #removed "gnomADg_MAF_afr", "gnomADg_AF_eas"
            out.write("###{}\t{}\t{}\t{}\t{}\t"
                  "{}\t{}\t{}\t{}\t{}\t"
                  "{}\t{}\t{}\t{}\t{}\t"
                  "{}\t{}\t{}\t{}\t{}\t"
                  "{}\t{}\t{}\t{}\t{}\t"
                  "{}\t{}\n".format('\t'.join(s), "Chrom", "Position", "Ref_allele", "Alt_allele", "rsID",
                                            "MAX_AF", "MAX_AF_POPS", "1000G", "gnomADg_MAF_global", "gnomADg_MAF_nfe",
                                            "gnomADg_MAF_nfe_nwe", "GERP(>3.6)",
                                            "phastCons(>0.71)",
                                            "phyloP(>0.8)", "CADD (>15)", "dbscSNV (>0.6)", "MaxEntScan (>1)",
                                            "Consequence_1",
                                            "Impact_1", "Gene_1", "Transcript_1", "HGVSg", "HGVSc",
                                            "HGVSp", "location", "Clinical_Sign_1", "Consequence_2..,etc"))

        else:
            s = ["Sample_ID", "isPositive"] if negative_samples_provided else ["Sample_ID"]
            out.write("###{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".
                      format('\t'.join(s), "Chrom", "Position", "Ref_allele", "Alt_allele", "rsID",
                            "MAX_AF", "MAX_AF_POPS", "1000G", "gnomADg_MAF_global", "gnomADg_MAF_nfe",
                            "gnomADg_MAF_nfe_nwe", "location"))
            for field in specific_fields:
                out.write("\t{}".format(field))
            out.write("\n")

        for sample, variantSet in sample_based_dict.items():

            s = [sample, "no"] if sample in negative_samples_provided else [sample, "yes"] if \
                negative_samples_provided else [sample]
            for variant in variantSet:
                out.write("{}\t{}\n".format('\t'.join(s), '\t'.join(str(v) for v in variant)))
    out.close()


def extractFields(vcf, output_file, negInd, fields):
    if negInd:
        with open(negInd,'r') as infile:
            negative_samples_provided = [line.strip() for line in infile]
    else:
        negative_samples_provided=[]

    vcf_data = VCF(vcf, gts012=True)
    samples = vcf_data.samples
    sample_based=defaultdict(list)
    hp = hgvs.parser.Parser()
    tools = []
    for field in vcf_data.header_iter():
        if field["HeaderType"] == "INFO" and field["ID"] == "ANN":
            tools = field["Description"].split("Format:")[1][:-1].strip().split("|")
            print(tools)

    for record in vcf_data:
        heterozygous = (record.gt_types == 1).nonzero()[0].tolist()
        hom_alt = (record.gt_types == 2).nonzero()[0].tolist()
        merged_individuals_with_variant=heterozygous+hom_alt
        if fields is not None:
            for ind in merged_individuals_with_variant:
                sample_based = get_annotations(record, tools, samples, sample_based, ind, hp, specific_fields=fields)

        else:
            for ind in merged_individuals_with_variant:
                try:
                    sample_based = get_annotations(record, tools, samples, sample_based, ind, hp)
                except KeyError:
                    print("Record {},{},{},{} does not have ANN field.".format(record.CHROM, str(record.POS),
                                                                               record.REF, record.ALT[0]))

        write_output(sample_based, output_file, negative_samples_provided, fields)


def main():
    parser = argparse.ArgumentParser(description='Script that produces sample based excel readable file for annotated variants.')
    parser.add_argument(dest='vcf', help='Path to the vcf. This file must have been split and normalized before.')
    parser.add_argument(dest='output_file', help='Path to the output tab separated file.')
    parser.add_argument('-n','--negativeIndividuals', metavar="neg_indivuals", help='File with negative sample names for a given phenotype. Will add column isPositive to the output.')
    parser.add_argument('-f','--fields', nargs='+', help='Specific fields to include')
    args = parser.parse_args()

    extractFields(args.vcf, args.output_file, args.negativeIndividuals, args.fields)


if __name__ == "__main__":
    main()
