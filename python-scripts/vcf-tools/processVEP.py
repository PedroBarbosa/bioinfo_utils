import numpy as np
import vcf
import argparse
import logging
import sys
import operator
import os
from collections import OrderedDict, defaultdict
logging.basicConfig(stream=sys.stdout, level=logging.INFO,  format='%(asctime)s %(message)s')


noTranscriptMatch = 0
transcriptMatch = 0
multipleTranscriptMatch = 0
totalVariants=0
passedVariants=0
partialPassed=0
allFailed=0
siftcounts=0
def updateTotalVariants():
    global totalVariants
    totalVariants +=1
    return totalVariants

def updateVariantsPassingFilters():
    global passedVariants
    passedVariants +=1
    return passedVariants

def updateVariantsPartiallyPassing():
    global partialPassed
    partialPassed +=1
    return partialPassed

def updateVariantsFailing():
    global allFailed
    allFailed +=1
    return allFailed
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

def checkValidFiltersSyntax(filtersFile,anno_fields,attributes=None,operation=False):

    if operation:
        ops=["greater","lower","equal","contains"]
        operations = [x.rstrip().split('\t')[2] for x in open(filtersFile).readlines()]
        op_values = [x.rstrip().split('\t')[1] for x in open(filtersFile).readlines()]
        #operations = np.genfromtxt(filtersFile, dtype=str, delimiter='\t', usecols=(2))
        #operations = np.atleast_1d(operations)
        for v in operations:
            for op in v.rstrip().split("_"):
                if op not in ops:
                    logging.error("ERROR. Invalid operation found in 3rd column of filter file: [{}]".format(op))
                    exit(1)

        for v in op_values:
            if "_" in v:
                for val in v.split("_"):
                    try:
                        float(val)
                        logging.error("ERROR. Multiple values (splitted by '_') per filter on the 2nd column ({}) are only available when comparing strings (e.g contains operation).".format(val))
                        exit(1)
                    except ValueError:
                        continue
    elif attributes:
        for attr in attributes[1]:

            if attr not in attributes[0]:
                if attr not in anno_fields:
                    print(attr)
                    logging.error("ERROR: {} filter atrribute is not present in VCF record. Please set valid filters in the INFO field names in the filters file".format(attr))
                    exit(1)
                else:
                    logging.info("{} attribute specified is located within the ANNO field.".format(attr))

def filtersToDict(filtersFile):
    #attribute=np.genfromtxt(filtersFile, dtype=str, usecols=(0))
    #operation = np.genfromtxt(filtersFile, dtype=str, delimiter='\t', usecols=range(1, 3))
    attribute = [x.split('\t')[0] for x in open(filtersFile).readlines()]
    operation=[x.rstrip().split('\t')[1:3] for x in open(filtersFile).readlines()]
    return OrderedDict(zip(attribute, operation))

def applyFilterWithinANNO(vcfrecord,attr,filters,ops,anno_fields,noneDiscard,justFirstConsequence,filter_result):
    passed_detected=False
    dict_multiple_cons=defaultdict(list)
    for consequences in vcfrecord.INFO["ANN"]:
        cons_fiels=consequences.split("|")

        #for i in range(0,len(anno_fields)):
        #    print(anno_fields[i],cons_fiels[i])
        if passed_detected == False:
            if not cons_fiels[anno_fields.index(attr)]:
                if noneDiscard == False:
                    filter_result.append((attr, "PASS_asNA"))
                    passed_detected=True
                else:
                    dict_multiple_cons[attr].append("FAIL_asNA")
                    #filter_result.append((attr, "FAIL_asNA"))
            else:
                try:
                    a=float(cons_fiels[anno_fields.index(attr)])
                    b=float(filters[0])
                    if ops[filters[1]](a,b):
                        filter_result.append((attr, "PASS"))
                        passed_detected=True
                    else:
                        dict_multiple_cons[attr].append("FAIL_asNA")
                        #filter_result.append((attr, "FAIL"))
                except ValueError:
                    try:
                        inlength=len(filter_result)
                        for value in filters[0].split("_"):
                            if ops[filters[1]](cons_fiels[anno_fields.index(attr)], value):
                                filter_result.append((attr, "PASS"))
                                passed_detected=True
                                break #substring found, no need to search more

                        if len(filter_result) == inlength: #no match for multiple filters, thus fail
                            dict_multiple_cons[attr].append("FAIL")
                            #filter_result.append((attr, "FAIL"))

                    except TypeError:
                        logging.error("Error. {} filter contains invalid operations ({}) given value obtained in the VCF record ({}) the value ({}) provided. Perhaps you are comparing strings? If so, use contains operation "
                                      "operation?".format(attr,filters[1],cons_fiels[anno_fields.index(attr)],filters[0]))
                        exit(1)
        if justFirstConsequence:
            break
    if not filter_result:
        filter_result.append((attr,dict_multiple_cons[attr][0]))

    return filter_result

def applyFilter(vcfrecord,filterDict,anno_fields,noneDiscard,permissive,justFirstConsequence,outfailed):
    ops = {"equal_lower": operator.le, "lower_equal" : operator.le, "greater_equal" : operator.ge, "equal_greater" : operator.ge
           , "lower" : operator.lt, "greater" : operator.gt, "equal" : operator.eq, "contains" : operator.contains}

    filter_result=[]
    for attr,filters in filterDict.items():
        if not attr in vcfrecord.INFO:
            filter_result=applyFilterWithinANNO(vcfrecord,attr,filters,ops,anno_fields,noneDiscard,justFirstConsequence,filter_result)
        else:
            if vcfrecord.INFO[attr] is None:
                if noneDiscard == False:
                    filter_result.append((attr,"PASS_asNA"))
                else:
                    filter_result.append((attr, "FAIL_asNA"))
            elif isinstance(vcfrecord.INFO[attr], list) and vcfrecord.INFO[attr][0] is None:
                if noneDiscard == False:
                    filter_result.append((attr, "PASS_asNA"))
                else:
                    filter_result.append((attr, "FAIL_asNA"))
            elif isinstance(vcfrecord.INFO[attr], list) and ops[filters[1]](float(vcfrecord.INFO[attr][0]), float(filters[0])):
                filter_result.append((attr, "PASS"))
            elif isinstance(vcfrecord.INFO[attr], float) and ops[filters[1]](vcfrecord.INFO[attr], float(filters[0])):
                filter_result.append((attr, "PASS"))
            else:
                filter_result.append((attr, "FAIL"))

    with open(outfailed + '_filteringOutput.tsv','a') as filtout:
        filtout.write(str(vcfrecord) + "\t" + '\t'.join([x[1] for x in filter_result]) + "\n")
    filtout.close()

    if all("PASS" in x[1] for x in filter_result):
        updateVariantsPassingFilters()
        return vcfrecord
    elif any("PASS" in x[1] for x in filter_result):
        updateVariantsPartiallyPassing()
        if permissive == True:
            return vcfrecord
    elif all("FAIL" in x[1] for x in filter_result):
        updateVariantsFailing()

def vcfreader(invcf,transcriptIDs,filterConfig,bedLotation,outbasename,noneValues,permissive,just1stConsequence):
    vcf_reader = vcf.Reader(filename=invcf)
    vcf_writer = vcf.Writer(open(outbasename + '_filtered.vcf', 'w'), vcf_reader)
    checkFileExists(outbasename + "_failedFilter.txt")
    anno_fields=[]
    filters=False
    transcripts=False
    bedLocation=False

    if not "ANN" in vcf_reader.infos.keys():
        logging.error("VCF record doesn't seem to have the required ANN field to process transcript consequences or filters usually found within such field.")
        exit(1)
    else:
        anno_fields = vcf_reader.infos["ANN"][3].split(":")[1].lstrip().split("|")
        print(anno_fields)
    if filterConfig:
        filters=True
        checkValidFiltersSyntax(filterConfig,[],operation=True)
        filters_dict=filtersToDict(filterConfig)
        checkValidFiltersSyntax(filters_dict, anno_fields,attributes=(list(vcf_reader.infos.keys()),list(filters_dict.keys())))
        with open(outbasename + '_filteringOutput.tsv', 'w') as failed:
            failed.write('#Variant_record' + "\t" + '\t'.join(list(filters_dict.keys())) + "\n")
        failed.close()
    if transcriptIDs:
        transcripts=True
        transcripts_list=transcriptIDsToList(transcriptIDs)

    for record in vcf_reader:
        updateTotalVariants()
        if not "*" in record.ALT:
            if filters:
                record=applyFilter(record,filters_dict,anno_fields,noneValues,permissive,just1stConsequence,outbasename)

            if record and transcripts :
                record=filterVEPtranscripts(record,transcripts_list,anno_fields,just1stConsequence)

            if record:
                vcf_writer.write_record(record)
        else:
            with open(outbasename + '_failedFilter.txt', 'a') as failed:
                failed.write("{}\t{}\n".format(str(record), 'ALT allele asterisk'))
            failed.close()

    vcf_writer.close()
    logging.info("{}\t{}".format("Total number of variants:",totalVariants))
    if filters:
        logging.info("{}\t{}".format("Number of variants passing all the filters:",passedVariants))
        if len(filters_dict.keys()) > 1:
            logging.info("{}\t{}".format("Number of variants passing at least one the filters:", partialPassed))
        logging.info("{}\t{}".format("Number of variants failing the filters:", allFailed))

    if transcripts:
        logging.info("\n{}\t{}".format("Number of variants filtered with match to transcriptIDs:",transcriptMatch))
        logging.info("{}\t{}".format("Number of variants filtered with NO match to transcriptIDs:", noTranscriptMatch))
        logging.info("{}\t{}".format("Number of variants filtered with multiple matches to transcriptIDs:", multipleTranscriptMatch))



def main():
    parser = argparse.ArgumentParser(
        description='Script to process VCF vep output useful data from a genome coverage analysis performed with GAT4K4. Matplotlib is required')
    parser.add_argument(dest='vcf',  help='Path to the vcf')
    parser.add_argument('-f','--filter', help='Auxiliar tab delimited file with the INFO attributes to filter. E.g: "AF 0.01    greater". Valid operations for '
                                              'the 3rd col: [greater,lower,equal,contains]. "_" may be used to apply multiple operations.')
    parser.add_argument('-n','--none', action='store_true', help='Flag to disable None values (.) for specific filter to count as passing. E.g prediction tools with None assignments should not count as passing filters. '
                                                                 'Default: None values are treated as passing the filters (e.g 1000G mafs')
    parser.add_argument('-p','--permissive', action='store_true', help='Flag to allow the script to accept any variant that passes at least 1 filter, Default: A variant must comply with ALL the filters supplied in order'
                                                                      'to pass the filtering stage.')
    parser.add_argument('-c','--firstConsequence', action='store_true', help='Flag to apply filters on just the first consequence. Default: Applies to all. Particularly important when using ANN filters')
    parser.add_argument('-t','--transcriptIDs',help='File with ensembl transcript IDs to keep variant consequence within the ANN field.')
    parser.add_argument('-l','--location', help='Try to classify each variant as deep intronic/intronic/spliceSite/exonic based on an input bed file.')
    parser.add_argument("-o", "--output", required=True, help='Basename to write output files.')
    args = parser.parse_args()

    if not args.filter and not args.transcriptIDs and not args.deepIntronic:
        logging.error("ERROR: No operation was set. Please specificy '-f', '-t', '-l' or  a combination of them.")
        exit(1)

    vcfreader(args.vcf,args.transcriptIDs,args.filter,args.location,args.output,args.none,args.permissive,args.firstConsequence)
if __name__ == "__main__":
    main()

