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



def extractFields(vcf,negInd,fields):
    if negInd:
        with open(negInd,'r') as infile:
            negative = [line.strip() for line in infile]
    else:
        negative=[]
    indexes=defaultdict(list)
    vcf_data = VCF(vcf, gts012=True)
    samples = vcf_data.samples
    no_ref_genotypes=(1,2)
    sample_based=defaultdict(list)
    hp = hgvs.parser.Parser()
    for field in vcf_data.header_iter():
        if field["HeaderType"] == "INFO" and field["ID"] == "ANN":
            tools = field["Description"].split("Format:")[1][:-1].strip().split("|")
            print(tools)

    for record in vcf_data:
        heterozygous = (record.gt_types == 1).nonzero()[0].tolist()
        hom_alt = (record.gt_types == 2).nonzero()[0].tolist()
        merged=heterozygous+hom_alt
        if leng(fields) > 0:
            for i in merged:
                try:
                    var=[record.CHROM,str(record.POS), record.REF, record.ALT[0],record.INFO["ANN"].split(",")[0].split("|")[tools.index('Existing_variation')]
                    sample_based[samples[i]].append(tuple(var))
                except KeyError:
                    print("Record {},{},{},{} does not have ANN field.".format(record.CHROM,str(record.POS), record.REF, record.ALT[0]))

        else:
            for i in merged:

                try:
                    var=[record.CHROM,str(record.POS), record.REF, record.ALT[0],record.INFO["ANN"].split(",")[0].split("|")[tools.index('Existing_variation')],
                    record.INFO["ANN"].split(",")[0].split("|")[tools.index('MAX_AF')],
                    record.INFO["ANN"].split(",")[0].split("|")[tools.index('MAX_AF_POPS')],
                    record.INFO["ANN"].split(",")[0].split("|")[tools.index('AF')],
                    record.INFO["ANN"].split(",")[0].split("|")[tools.index('gnomADg_AF')],
                    record.INFO["ANN"].split(",")[0].split("|")[tools.index('gnomADg_AF_nfe')],
                    record.INFO["ANN"].split(",")[0].split("|")[tools.index('gnomADg_AF_nfe_nwe')],
                    #record.INFO["ANN"].split(",")[0].split("|")[tools.index('gnomADg_AF_afr')],
                    #record.INFO["ANN"].split(",")[0].split("|")[tools.index('gnomADg_AF_eas')],
                    record.INFO["ANN"].split(",")[0].split("|")[tools.index('GERP')],
                    record.INFO["ANN"].split(",")[0].split("|")[tools.index('phastCons')],
                    record.INFO["ANN"].split(",")[0].split("|")[tools.index('phyloP')],
                    record.INFO["ANN"].split(",")[0].split("|")[tools.index('CADD_PHRED')],
                    max(record.INFO["ANN"].split(",")[0].split("|")[tools.index('ada_score')],record.INFO["ANN"].split(",")[0].split("|")[tools.index('rf_score')]),
                    record.INFO["ANN"].split(",")[0].split("|")[tools.index('MaxEntScan_diff')]]

                    for cons in record.INFO["ANN"].split(","):
                        fields=cons.split("|")
                        var.extend((fields[tools.index('Consequence')],fields[tools.index('IMPACT')],fields[tools.index('SYMBOL')],
                                fields[tools.index('Feature')],
                                fields[tools.index('HGVSg')],
                                fields[tools.index('HGVSc')],
                                fields[tools.index('HGVSp')],
                                get_location(fields[tools.index('HGVSc')],hp),
                                fields[tools.index('CLIN_SIG')],
                                ))

                    sample_based[samples[i]].append(tuple(var))
                except KeyError:
                    print("Record {},{},{},{} does not have ANN field.".format(record.CHROM,str(record.POS), record.REF, record.ALT[0]))



    with open('excel.tsv','w') as out:
        s = ["Sample_ID","isPositive"] if negative else ["Sample_ID"]
        out.write("###{}\t{}\t{}\t{}\t{}\t"
                  "{}\t{}\t{}\t{}\t{}\t"
                  "{}\t{}\t{}\t{}\t{}\t"
                  "{}\t{}\t{}\t{}\t{}\t"
                  "{}\t{}\t{}\t{}\t{}\t"
                  "{}\t{}\t{}\t{}\n".format('\t'.join(s),"Chrom","Position","Ref_allele","Alt_allele","rsID",
                    "MAX_AF","MAX_AF_POPS","1000G","gnomeADg_MAF_global","gnomeADg_MAF_nfe",
                    "gnomeADg_MAF_nfe_nwe","gnomADg_MAF_afr", "gnomADg_AF_eas","GERP(>3.6)","phastCons(>0.71)",
                    "phyloP(>0.8)","CADD (>15)","dbscSNV (>0.6)","MaxEntScan (>1)","Consequence_1",
                    "Impact_1","Gene_1","Transcript_1","HGVSg","HGVSc",
                    "HGVSp","location","Clinical_Sign_1","Consequence_2..,etc"))
        
        for sample,variantSet in sample_based.items():

            s = [sample,"no"] if sample in negative else [sample,"yes"] if negative else [sample]
            for variant in variantSet:
                out.write("{}\t{}\n".format('\t'.join(s),'\t'.join(variant)))



def main():
    parser = argparse.ArgumentParser(description='Script to produces excel readable file for annotated variants.')
    parser.add_argument(dest='vcf', help='Path to the vcf. This file mast have been split and normalized before.')
    parser.add_argument('-n','--negativeIndividuals', help='File with negative sample names for a given phenotype. Will add column isPositive to the output.')
    parser.add_argument('-f','--fields', nargs='+', help='Additional fields to include')
    args = parser.parse_args()

    extractFields(args.vcf,args.negativeIndividuals,args.fields)
if __name__ == "__main__":
    main()
