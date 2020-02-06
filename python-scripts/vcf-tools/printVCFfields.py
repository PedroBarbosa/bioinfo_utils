import argparse
from cyvcf2 import VCF
from collections import defaultdict

def gets_standard_field(record, field):
    if field == "CHROM":
        return record.CHROM
    elif field == "POS":
        return str(record.POS)
    elif field == "ID":
        return str(record.ID)
    elif field == "REF":
        return record.REF
    elif field == "ALT":
        return ','.join(record.ALT)
    elif field == "QUAL":
        return str(record.QUAL)
    elif field == "FILTER":
        return record.FILTER
    elif field == "FORMAT":
        return ','.join(record.FORMAT)
    else:
        if record.INFO.get(field) == None:
            return "."
        else :
            return str(record.INFO.get(field))


def printFields(vcf,fields,printall):
    standard_fields=["CHROM","POS","ID","REF","ALT","QUAL","FILTER","FORMAT"]
    existing_info, final_final_existing_info=[],[]
    indexes=defaultdict(list)
    vcf_data = VCF(vcf, gts012=True)
    for field in vcf_data.header_iter():
        d = field.info()
        if d['HeaderType'] == "INFO":
            existing_info.append(d['ID'])
    print("#List of available INFO fields extracted from the header:\n#{}".format(existing_info))
    for field in vcf_data.header_iter():
        if field["HeaderType"] == "INFO" and field["ID"] == "ANN":
            tools = field["Description"].split("Format:")[1][:-1].strip().split("|")
            print("#List of available fields within ANN field:\n#{}".format(tools))
            final_existing_info = list(set(existing_info) - set(tools))
            for f in fields:
                try:
                    if not f in standard_fields and not f in final_existing_info:
                        indexes[f] = tools.index(f)
                except ValueError:
                    print("{} field not recognized".format(f))
                    exit(1)
        else:
            final_existing_info = existing_info
    
    print("#{}".format("\t".join(fields)))
    for record in vcf_data:
        d=defaultdict(list)
        outline=[]
        try:
            info_obj = record.INFO
        except AttributeError:
            standard_where_no_info=[]
            [standard_where_no_info.append(gets_standard_field(record,f)) for f in fields if f in standard_fields]
            if len(standard_where_no_info) > 0:
                print('\t'.join(standard_where_no_info))
            continue
        
        try:
            if printall:
                info_vep = info_obj.get("ANN").split(",")
            else:
                info_vep = info_obj.get("ANN").split(",")[0].split("|")
        except AttributeError:
            info_vep=[]
        
        for f in fields:
            if printall and not info_vep:
                print("No ANN field found in VCF, therefore --printAll argument is useless.\n")
                exit(1)
            elif printall:
                for i,block in enumerate(info_vep):
                    try:
                        d[i].append(gets_standard_field(record,f)) if f in standard_fields or f in final_existing_info else d[i].append(block.split("|")[indexes[f]])
                    except IndexError:
                        continue
            else:
                if f in standard_fields or f in final_existing_info:
                    outline.append(gets_standard_field(record,f))
                elif info_vep:
                    [outline.append(info_vep[indexes[f]]) if info_vep[indexes[f]] else outline.append('None')]
                else:
                    print("{} field not in VCF".format(f))
                    exit(1)
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

