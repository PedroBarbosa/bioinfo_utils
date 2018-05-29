import argparse
from cyvcf2 import VCF
from collections import defaultdict

def printFields(vcf,fields):
    indexes=defaultdict(list)
    vcf_data = VCF(vcf, gts012=True)
    for field in vcf_data.header_iter():
        if field["HeaderType"] == "INFO" and field["ID"] == "ANN":
            tools = field["Description"].split("Format:")[1][:-1].strip().split("|")
            print(tools)
            for f in fields:
                try:
                    indexes[f] = tools.index(f)
                except ValueError:
                    print("{} field not recognized".format(f))
                    exit(1)
    print("\t".join(indexes.keys()))
    for record in vcf_data:
        outline=[]
        try:
            info = record.INFO.get("ANN").split(",")[0].split("|")
        except AttributeError:
            continue
        for f in fields:
            [outline.append(info[indexes[f]]) if info[indexes[f]] else outline.append('None')]
        print('\t'.join(outline))

def main():
    parser = argparse.ArgumentParser(description='Script print specific fields from VCF to tab delimited')
    parser.add_argument(dest='vcf', help='Path to the vcf')
    parser.add_argument(dest='fields', nargs='+', help='Fields to parse')
    args = parser.parse_args()

    printFields(args.vcf,args.fields)
if __name__ == "__main__":
    main()
