import argparse
from cyvcf2 import VCF
from collections import defaultdict

def extractFields(vcf,fields):
    indexes=defaultdict(list)
    vcf_data = VCF(vcf, gts012=True)
    samples = vcf_data.samples
    no_ref_genotypes=(1,2)
    sample_based=defaultdict(list)

    for field in vcf_data.header_iter():
        if field["HeaderType"] == "INFO" and field["ID"] == "ANN":
            tools = field["Description"].split("Format:")[1][:-1].strip().split("|")
            print(tools)

    for record in vcf_data:
        heterozygous = (record.gt_types == 1).nonzero()[0].tolist()
        hom_alt = (record.gt_types == 2).nonzero()[0].tolist()
        merged=heterozygous+hom_alt
        for i in merged:
            var=[record.CHROM,str(record.POS), record.REF, record.ALT[0],record.INFO["ANN"].split(",")[0].split("|")[tools.index('Existing_variation')],
                 record.INFO["ANN"].split(",")[0].split("|")[tools.index('AF')],
                 record.INFO["ANN"].split(",")[0].split("|")[tools.index('gnomADg_AF_NFE')],
                 record.INFO["ANN"].split(",")[0].split("|")[tools.index('GERP')],
                 record.INFO["ANN"].split(",")[0].split("|")[tools.index('phastCons')],
                 record.INFO["ANN"].split(",")[0].split("|")[tools.index('phyloP')]]
            for cons in record.INFO["ANN"].split(","):
                fields=cons.split("|")
                var.extend((fields[tools.index('Consequence')],fields[tools.index('IMPACT')],fields[tools.index('SYMBOL')],
                            fields[tools.index('Feature')],
                            fields[tools.index('HGVSg')],
                            fields[tools.index('HGVSc')],
                            fields[tools.index('HGVSp')],
                            fields[tools.index('CLIN_SIG')],
                            ))

            sample_based[samples[i]].append(tuple(var))


    with open('excel.tsv','w') as out:
        out.write("###{}\t{}\t{}\t{}\t{}\t"
                  "{}\t{}\t{}\t{}\t{}\t"
                  "{}\t{}\t{}\t{}\t{}\t"
                  "{}\t{}\t{}\t{}\t{}\n".format("Sample_ID","Chrom","Position","Ref_allele","Alt_allele",
                                                "rsID", "1000G","gnomeAD_MAF","GERP(>3.6)","phastCons(>0.71)",
                                                "phyloP(>0.8)","Consequence_1","Impact_1", "Gene_1","Transcript_1",
                                                "HGVSg","HGVSc","HGVSp","Clinical_Sign_1","Consequence_2..,etc"))
        for sample,variantSet in sample_based.items():
            for variant in variantSet:
                out.write("{}\t{}\n".format(sample,'\t'.join(variant)))



def main():
    parser = argparse.ArgumentParser(description='Script to produces excel readable file for annotated variants.')
    parser.add_argument(dest='vcf', help='Path to the vcf. This file mast have been split and normalized before.')
    parser.add_argument('-f','--fields', nargs='+', help='Additional fields to include')
    args = parser.parse_args()

    extractFields(args.vcf,args.fields)
if __name__ == "__main__":
    main()
