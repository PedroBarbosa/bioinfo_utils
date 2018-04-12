import numpy as np
import vcf
import argparse
import logging
import sys
import operator
import os
logging.basicConfig(stream=sys.stdout, level=logging.INFO,  format='%(asctime)s %(message)s')


noTranscriptMatch = 0
transcriptMatch = 0
multipleTranscriptMatch = 0
totalVariants=0
passedVariants=0

def updateTotalVariants():
    global totalVariants
    totalVariants +=1
    return totalVariants

def updateVariantsPassingFilters():
    global passedVariants
    passedVariants +=1
    return passedVariants

def updateMultipleTranscriptMatch():
    global multipleTranscriptMatch
    multipleTranscriptMatch+=1
    return multipleTranscriptMatch

def updateNoTranscriptMatch():
    global noTranscriptMatch
    noTranscriptMatch+=1
    return noTranscriptMatch

def updateTranscriptMatch():
    global transcriptMatch
    transcriptMatch+=1
    return transcriptMatch

def checkFileExists(file):
    try:
        os.remove(file)
    except OSError:
        pass
def transcriptIDsToList(transcripts):
    with open(transcripts,"r") as infile:
        return [line.split(".")[0] for line in infile]

def filterVEPtranscripts(vcfrecord,transcripts,anno_fields):






    #print(anno_fields)
    if not "ANN" in vcfrecord.INFO.keys():
        logging.error("VCF record doesn't seem to have the required ANN field to process transcript consequences.")
        print(vcfrecord.ALT)
        exit(1)
    else:
        consequences=vcfrecord.INFO["ANN"]
        new_ann=[]
        new_variant=True
        for cons in consequences:
            cons_fields=cons.split("|")
            if cons_fields[anno_fields.index("Feature_type")] == "Transcript" and cons_fields[anno_fields.index("Feature")] in transcripts:
                new_ann.append("|".join(cons_fields))
                if new_variant:
                    updateTranscriptMatch()
                    new_variant=False

        if len(new_ann) == 0:
            updateNoTranscriptMatch()
            return vcfrecord
        elif len(new_ann) == 1:
            vcfrecord.INFO["ANN"] = new_ann
            return vcfrecord
        else:
            updateMultipleTranscriptMatch()
            vcfrecord.INFO["ANN"] = new_ann
            return vcfrecord




def checkValidFiltersSyntax(filtersFile,attributes=None,operation=False):
    if operation:
        ops=["greater","lower","equal"]
        operations = np.genfromtxt(filtersFile, dtype=str, delimiter='\t', usecols=(2))
        for v in operations:
            for op in v.rstrip().split("_"):
                if op not in ops:
                    logging.error("ERROR. Invalid operation found in 3rd column of filter file: [{}]".format(op))
                    exit(1)
    elif attributes:
        for attr in attributes[1]:
            if attr not in attributes[0]:
                logging.error("ERROR: {} filter atrribute is not present in VCF record. Please set valid filters in the INFO field names in the filters file".format(attr))
                exit(1)

def filtersToDict(filtersFile):
    attribute=np.genfromtxt(filtersFile, dtype=str, usecols=(0))
    operation = np.genfromtxt(filtersFile, dtype=str, delimiter='\t', usecols=range(1, 3))
    return dict(zip(attribute, operation.tolist()))

def applyFilter(vcfrecord,filterDict,outfailed):
    checkValidFiltersSyntax(filterDict,attributes=(list(vcfrecord.INFO.keys()),list(filterDict.keys())))
    ops = {"equal_lower": operator.le, "lower_equal" : operator.le, "greater_equal" : operator.ge, "equal_greater" : operator.ge
           , "lower" : operator.lt, "greater" : operator.gt, "equal" : operator.eq}

    filter_result=[]
    #for field in vcfrecord.INFO:
        #if field in filterDict.keys():

    for attr,filters in filterDict.items():
        #print(vcfrecord.INFO)
        #print(attr,vcfrecord.INFO[attr])
        #if vcfrecord.POS == 116275093:
        #    debug=1
        #at=vcfrecord.INFO[attr]
        if vcfrecord.INFO[attr] is None:
            filter_result.append((attr,"PASS"))
        elif isinstance(vcfrecord.INFO[attr], list) and vcfrecord.INFO[attr][0] is None:
            filter_result.append((attr, "PASS"))
        elif isinstance(vcfrecord.INFO[attr], list) and ops[filters[1]](float(vcfrecord.INFO[attr][0]), float(filters[0])):
            filter_result.append((attr, "PASS"))
        elif isinstance(vcfrecord.INFO[attr], float) and ops[filters[1]](vcfrecord.INFO[attr], float(filters[0])):
            filter_result.append((attr, "PASS"))
        else:
            filter_result.append((attr, "FAIL"))

    if all(x[1] == "PASS" for x in filter_result):
        updateVariantsPassingFilters()
        return vcfrecord
    else:
        #print(vcfrecord,filter_result)
        with open(outfailed + '_failedFilter.txt','a') as failed:
            failed.write("{}\t{}\n".format(str(vcfrecord),','.join([x[0] for x in filter_result if x[1] != "PASS"])))
        failed.close()



def vcfreader(invcf,transcriptIDs,filterConfig,bedLotation,outbasename):
    vcf_reader = vcf.Reader(filename=invcf)
    vcf_writer = vcf.Writer(open(outbasename + '_filtered.vcf', 'w'), vcf_reader)
    checkFileExists(outbasename + "_failedFilter.txt")
    filters=False
    transcripts=False
    bedLocation=False
    if filterConfig:
        filters=True
        checkValidFiltersSyntax(filterConfig,operation=True)
        filters_dict=filtersToDict(filterConfig)

    if transcriptIDs:
        transcripts=True
        transcripts_list=transcriptIDsToList(transcriptIDs)
    anno_fields=vcf_reader.infos["ANN"][3].split(":")[1].lstrip().split("|")
    for record in vcf_reader:
        updateTotalVariants()
        if not "*" in record.ALT:
            if filters:
                record=applyFilter(record,filters_dict,outbasename)

            if record and transcripts :
                record=filterVEPtranscripts(record,transcripts_list,anno_fields)

            if record:
                vcf_writer.write_record(record)
        else:
            with open(outbasename + '_failedFilter.txt', 'a') as failed:
                failed.write("{}\t{}\n".format(str(record), 'ALT allele asterisk'))
            failed.close()

    vcf_writer.close()

    logging.info("{}\t{}".format("Total number of variants:",totalVariants))
    if filters:
        logging.info("{}\t{}\n".format("Number of variants passing the filters:",passedVariants))
    if transcripts:
        logging.info("{}\t{}".format("Number of variants filtered with match to transcriptIDs:",transcriptMatch))
        logging.info("{}\t{}".format("Number of variants filtered with NO match to transcriptIDs:", noTranscriptMatch))
        logging.info("{}\t{}".format("Number of variants filtered with multiple matches to transcriptIDs:", multipleTranscriptMatch))


def main():
    parser = argparse.ArgumentParser(
        description='Script to process VCF vep output useful data from a genome coverage analysis performed with GAT4K4. Matplotlib is required')
    parser.add_argument(dest='vcf',  help='Path to the vcf')
    parser.add_argument('-f','--filter', help='Auxiliar tab delimited file with the INFO attributes to filter. E.g: "AF 0.01    greater". Valid operations for '
                                              'the 3rd col: [greater,lower,equal]. "_" may be used to apply multiple operations.')
    parser.add_argument('-t','--transcriptIDs',help='File with ensembl transcript IDs to keep variant consequence within the ANN field.')
    parser.add_argument('-l','--location', help='Try to classify each variant as deep intronic/intronic/spliceSite/exonic based on an input bed file.')
    parser.add_argument("-o", "--output", required=True, help='Basename to write output files.')
    args = parser.parse_args()

    if not args.filter and not args.transcriptIDs and not args.deepIntronic:
        logging.error("ERROR: No operation was set. Please specificy '-f', '-t', '-l' or  a combination of them.")
        exit(1)

    #declareGlobal()
    vcfreader(args.vcf,args.transcriptIDs,args.filter,args.location,args.output)
if __name__ == "__main__":
    main()

