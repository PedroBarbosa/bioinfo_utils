import argparse
from cyvcf2 import VCF, Writer

def addField(vcf,name,value,description,outfile):
    vcf_data = VCF(vcf, gts012=True)
    vcf_data.add_info_to_header({'ID':  name, 'Description': description, 'Type':'String', 'Number': '1'})
    w = Writer(outfile, vcf_data)
    for record in vcf_data:
        record.INFO[name] = value
        w.write_record(record)
    w.close();vcf_data.close()

def main():
    parser = argparse.ArgumentParser(description='Script to add simple annotations to the INFO field of a VCF')
    parser.add_argument(dest='vcf', help='Path to the vcf')
    parser.add_argument(dest='field_name', help='Name of the field to add')
    parser.add_argument(dest='field_value',help='Value of the field to add')
    parser.add_argument(dest='field_description', help='Info to add in the header')
    parser.add_argument(dest='outfile', help='Name of the output file')
    args = parser.parse_args()

    addField(args.vcf,args.field_name, args.field_value, args.field_description, args.outfile)
if __name__ == "__main__":
    main()