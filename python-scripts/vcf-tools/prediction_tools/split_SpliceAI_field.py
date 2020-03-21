import argparse
from cyvcf2 import VCF, Writer
import os


def set_new_fields(record, spliceAI_pred):

    gene, ds_ag, ds_al, ds_dg, ds_dl = "", "", "", "", ""
    predictions = spliceAI_pred.split(",")
    for i, prediction in enumerate(predictions):
        pred = prediction.split("|")
        if i != 0:
            gene = gene + "," + pred[1]
            ds_ag = ds_ag + "," + pred[2]
            ds_al = ds_al + "," + pred[3]
            ds_dg = ds_dg + "," + pred[4]
            ds_dl = ds_dl + "," + pred[5]
        else:
            gene = pred[1]
            ds_ag = pred[2]
            ds_al = pred[3]
            ds_dg = pred[4]
            ds_dl = pred[5]

    record.INFO['Gene_SpliceAI'] = gene
    record.INFO['DS_AG'] = ds_ag
    record.INFO['DS_AL'] = ds_al
    record.INFO['DS_DG'] = ds_dg
    record.INFO['DS_DL'] = ds_dl
    return record


def process_vcf(vcf):

    if not os.path.isfile(vcf):
        print("{} is not a valid file. Exiting.".format(vcf))
        exit(1)
    vcf_data = VCF(vcf, gts012=True)
    out = str(vcf).replace(".vcf", "_spliceAI_processed.vcf")

    vcf_data.add_info_to_header({'ID': 'Gene_SpliceAI',
                                 'Description': 'Gene for which spliceAI gave the prediction.',
                                 'Type': 'String', 'Number': '.'})
    vcf_data.add_info_to_header({'ID': 'DS_AG',
                                 'Description': 'SpliceAI score for an acceptor gain.',
                                 'Type': 'String', 'Number': '.'})
    vcf_data.add_info_to_header({'ID': 'DS_AL',
                                 'Description': 'SpliceAI score for an acceptor lost.',
                                 'Type': 'String', 'Number': '.'})
    vcf_data.add_info_to_header({'ID': 'DS_DG',
                                 'Description': 'SpliceAI score for a donor gain.',
                                 'Type': 'String', 'Number': '.'})
    vcf_data.add_info_to_header({'ID': 'DS_DL',
                                 'Description': 'SpliceAI score for a donor lost.',
                                 'Type': 'String', 'Number': '.'})
    w = Writer(out, vcf_data)
    for record in vcf_data:
        snvs = record.INFO.get('SpliceAI')
        indels = record.INFO.get('SpliceAI_ind')
        if snvs:
            record = set_new_fields(record, snvs)
        elif indels:
            record = set_new_fields(record, indels)
        w.write_record(record)

    vcf_data.close()
    w.close()
def main():
    parser = argparse.ArgumentParser(description='Script to split SpliceAI field into separate scores. Useful to'
                                                 'later prioritize variants based on specific thresholds.')
    parser.add_argument(dest='vcf', help='Path to the vcf')
    args = parser.parse_args()
    process_vcf(args.vcf)


if __name__ == "__main__":
    main()

