import argparse
from cyvcf2 import VCF
from collections import defaultdict

def printFields(vcf,fields,printall):
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
            if printall:
                info = record.INFO.get("ANN").split(",")
            else:
                info = record.INFO.get("ANN").split(",")[0].split("|")

        except AttributeError:
            continue
        d=defaultdict(list)
        for f in fields:
            try:
                if printall:
                    for i,block in enumerate(info):
                        d[i].append(block.split("|")[indexes[f]])

                else:
                    [outline.append(info[indexes[f]]) if info[indexes[f]] else outline.append('None')]
            except AttributeError:
                continue

        if len(d) > 0:
            for block in d.keys():
                print('\t'.join(d[block]))
            print('\n')
        else:
            print('\t'.join(outline))

def main():
    parser = argparse.ArgumentParser(description='Script print specific fields from VCF to tab delimited')
    parser.add_argument(dest='vcf', help='Path to the vcf')
    parser.add_argument(dest='fields', nargs='+', help='Fields to parse')
    parser.add_argument('-a', '--allConsequences', action='store_true',
                        help='Flag to print the values for all consequences blocks. Default: print the first.')
    args = parser.parse_args()

    printFields(args.vcf,args.fields, args.allConsequences)
if __name__ == "__main__":
    main()
