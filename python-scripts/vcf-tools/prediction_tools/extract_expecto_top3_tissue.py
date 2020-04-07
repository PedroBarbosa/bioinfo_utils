import argparse
import os
from cyvcf2 import VCF, Writer
import pysam


def checkIndex(expectofile):
    if not os.path.isfile(expectofile + '.tbi'):
        print('Index file is required. Make sure that {} file exists in the same path'.format(expectofile + '.tbi'))
        exit(1)


def processVCF(invcf,out,expecto):
    vcf_data = VCF(invcf, gts012=True)
    tbx_expecto = pysam.TabixFile(expecto)
    vcf_data.add_info_to_header({'ID': 'ExPecto', 'Description': 'A deep neural network aimed to predict the efffect of'
                                                                 ' variants on gene expression', 'Type': 'String',
                                 'Number': '.'})
    w = Writer(out, vcf_data)
    for record in vcf_data:
        try:
            for row in tbx_expecto.fetch(record.CHROM, record.start, record.end):
                if int(row.split()[1]) == record.POS and row.split()[2] == record.REF and row.split()[3] == record.ALT[0]:
                    a=1
                    #record.INFO["ExPecto"] = round(float(row.split()[4]), 2)
                    #break
                else:
                    record.INFO["DANN"] = "."
        except ValueError:
            record.INFO["ExPecto"] = "."

        w.write_record(record)


def main():
    parser = argparse.ArgumentParser(description='Script to aid on the task of add ExPecto scores to VCF files.')
    parser.add_argument(dest='vcf', default="-", help='Path to the vcf (uncompressed)')
    parser.add_argument(dest='out', help='Path to the VCF output file')
    parser.add_argument('--path', '--expectoPath', default="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/expecto",
                        help="Path where ExPecto files are located")
    parser.add_argument("-n", "--number", default=3, help="Top n number of tissues to report highest expression effects."
                                                          "Default:3")
    parser.add_argument("-d", "--directionality", action="store_true", help="If set, directionality scores will be included")
    args = parser.parse_args()

    expecto_filename = "expecto_1kb_TSS_v0.3.txt.gz"
    expecto_fullpath = os.path.join(args.expectoPath, expecto_filename)
    if os.path.isdir(args.expectoPath) and os.path.isfile(expecto_fullpath):
        checkIndex(expecto_fullpath)
        processVCF(args.vcf, args.out, expecto_fullpath, args.directionality)

    else:
        print("Does the ExPecto directory provided exists? Also, make sure the expecto scores (expecto_1kb_TSS_v0.3.txt.gz)"
              "are stored there.")
        exit(1)


if __name__ == "__main__":
    main()