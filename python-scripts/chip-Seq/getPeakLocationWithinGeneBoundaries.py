import argparse
import sys
import logging
import operator
from collections import defaultdict
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')


def readGff(inputgff):
    dict = defaultdict(list)

    with open(inputgff[0], 'r') as file:
        for line in file:
            line.rstrip()

            if not line.startswith("#"):
                feature = line.split("\t")
                if feature[2] == "gene":
                    gene_id = feature[8].rstrip()
                    dict[gene_id] = []
                elif feature[2] != "transcript":
                    feature_tuple = (feature[2], feature[3], feature[4], feature[6])
                    dict[gene_id].append(feature_tuple)
    file.close()
    return dict


def processPeaksFile(peakFile,gffDict, outputFile):

    with open(peakFile[0], 'r') as peakFile:
        with open(outputFile, 'w') as outFile:
            for line in peakFile:
                line =line.rstrip()
                if line.startswith("#"):
                    outFile.write(line + '\n')
                elif not "total_within" in line:
                    outFile.write(line + '\n')
                else:
                    fields = line.split("\t")
                    tss,tts,start_codon,stop_codon=False,False,False,False
                    i = fields.index("total_within")
                    gene_id = fields[i-2]
                    start_peak = fields[2]
                    end_peak = fields[3]
                    list_features = gffDict[gene_id]
                    result=[]
                    outside=False
                    outside_comment=""
                    max_coord=0

                    min_coord = min(list_features, key = lambda t: int(t[1]))[1]

                    for feature in list_features:
                        if int(feature[2]) > max_coord:
                            max_coord = int(feature[2])

                        if feature[1] == feature[2]: #if tts or tss
                            if feature[3] == "+":
                                if int(start_peak) <= int(feature[1]) and int(end_peak) >= int(feature[1]) and feature[0] == "tss":
                                    tss=True
                                elif feature[0] == "tts" and stop_codon == True:
                                    tts=True
                                    result.append("3'prime UTR")
                                elif feature[0] == "tts" and stop_codon == False:
                                    tts=True
                                    #result.append("3'prime UTR - without stop_codon")
                                    a = "cant' say nothin"

                            elif feature[3] == "-":
                                if feature[1] == "tts":
                                    tts=True
                                elif feature[0] == "tss" and start_codon == True:
                                    result.append("5'prime UTR")
                                    tss=True
                                elif feature[0] == "tss" and start_codon == False:
                                    tss=True
                                    #result.append("5'prime UTR - without start codon")
                                    a= "can't say nothing"



                        elif int(start_peak) >= int(feature[1]) and int(end_peak) <= int(feature[2]): #totally included in feature
                            result.append(feature[0])

                        elif int(start_peak) < int(feature[1]) and int(end_peak) > int(feature[2]) and not "codon" in feature[0] : #totally spans a feature, that is not stop/start codon and tts/tts
                            result.append(feature[0])

                        elif int(start_peak) > int(feature[1]) and int(start_peak) < int(feature[2]) and not "codon" in feature[0]: #partially included in the end of feature
                            result.append(feature[0])

                        elif int(start_peak) < int(feature[1]) and int(end_peak) > int(feature[1]) and not "codon" in feature[0]: #partially inclided in the beginning of feature
                            result.append(feature[0])

                        if feature[3] == "+": #tss appears before and tss appears after
                            if tss and "start_codon" in feature[0]:

                                if int(start_peak) < int(feature[1]) and int(end_peak) > int(feature[2]): # if totally included

                                    result.append("5'prime UTR")

                                elif int(start_peak) < int(feature[1]) and int(end_peak) < int(feature[1]): #if totally upstream
                                    result.append("5'prime UTR")

                            elif "start_codon" in feature[0]:
                                if int(start_peak) < int(feature[1]) and int(end_peak) > int(feature[2]): # if totally included
                                    result.append("5'prime UTR - without tss")
                                elif int(start_peak) < int(feature[1]) and int(end_peak) < int(feature[1]): #if totally upstream
                                    result.append("5'prime UTR - without tss")



                            elif "stop_codon" in feature[0]:
                                if int(start_peak) <= int(feature[1]) and int(end_peak) >= int(feature[2]):
                                    stop_codon=True



                        elif feature[3] == "-": #tts appears before
                            if tts and "stop_codon" in feature[0]:
                                if int(start_peak) < int(feature[1]) and int(end_peak) > int(feature[2]):
                                    result.append("3'prime UTR")
                                elif int(start_peak) < int(feature[1]) and int(end_peak) < int(feature[1]):
                                    result.append("3'prime UTR")
                                    #result.append("3'prime UTR - without tts")

                            elif "stop_codon" in feature[0]:
                                if int(start_peak) < int(feature[1]) and int(end_peak) > int(feature[2]): # if totally included
                                    result.append("3'prime UTR - without tss")
                                elif int(start_peak) < int(feature[1]) and int(end_peak) < int(feature[1]): #if totally upstream
                                    result.append("3'prime UTR - without tss")


                            elif "start_codon" in feature[0]:
                                if int(start_peak) <= int(feature[1]) and int(end_peak) >= int(feature[2]):
                                    start_codon = True




                        if int(end_peak) > max_coord and '+' in feature:
                            outside=True
                            outside_comment = "partial outside at 3' end (gene forward)"
                        elif int(start_peak) < int(min_coord) and '-' in feature:
                            outside= True
                            outside_comment = "partial outside at 3' end (gene reverse)"
                        else:
                            outside=False

                    if outside == True:
                        result.append(outside_comment)
                        logging.info("Peak %s, which has been marked totally within the gene %s, is actually partially outside of the gene at the downstream level" % (fields[0],gene_id))

                    if not result:
                        result =  "['No feature predicted']"
                        #logging.info("There is no feature detected in GFF within the gene %s, for which the peak %s was found to be totally within" % (gene_id,fields[0]))

                    fields[i] = fields[i] + ' ' + str(result)
                    outFile.write('\t'.join(fields) + '\n')
                    #print(fields)



def main():
    parser = argparse.ArgumentParser(description='Script to analyse peak location when it is located within the gene boundaries.')
    parser.add_argument(dest='peak_file', metavar='peakFile', nargs=1,
                        help='Excel file generated by the foldEnrichmentAnalysisFromMacs2 script with information peak location regarding the closest gene forward and reverse.')
    parser.add_argument(dest='gff', metavar='gffFile', nargs=1,
                        help="Genome GFF3 file to look deeper in the annotations.")
    parser.add_argument("-o", metavar='--outputFile', required=True, help="Output file to write the new excel.")
    args = parser.parse_args()

    logging.info("Processing GFF file..")
    gff_dict = readGff(args.gff)
    processPeaksFile(args.peak_file,gff_dict, args.o)

if __name__ == "__main__":
    main()
