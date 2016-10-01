import argparse
import sys
import logging
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

    ntss = 0
    ntts = 0
    ngenes = 0
    with open(peakFile[0], 'r') as peakFile:
        with open(outputFile, 'w') as outFile:
            for line in peakFile:
                line =line.rstrip()
                if line.startswith("#"):
                    outFile.write(line + '\n')
                elif not "total_within" in line:
                    outFile.write(line + '\n')
                elif "total_within" in line:
                    fields = line.split("\t")
                    i = fields.index("total_within")
                    gene_id = fields[i-2]
                    summit_peak = int(fields[5])
                    list_features = gffDict[gene_id]
                    result=[]
                    tts,tss,utr5,utr3, other_processed,start_cod, end_cod = False,False,False,False,False,False,False

                    ngenes +=1
                    min_coord = int(min(list_features, key = lambda t: int(t[1]))[1])
                    max_coord = int(max(list_features, key = lambda t: int(t[1]))[2])



                    for feature in list_features:
                        ######################################## FORWARD GENES ####################################333
                        if feature[3] == "+":
                            if int(feature[1]) == int(feature[2]) and int(feature[1]) == min_coord and tss== False and other_processed == False:
                                if feature[0] == "tss":
                                    tss = True
                                    ntss += 1
                            elif int(feature[1]) == int(feature[2]) and int(feature[1]) == max_coord and other_processed == True:
                                if feature[0] == "tts":
                                    tts = True
                                    ntts +=1
                            elif int(feature[1]) == int(feature[2]) and int(feature[1]) != min_coord and int(feature[1]) != max_coord :
                                logging.info("Warning. Gene id %s does not have tts or tss starting either at the end or beginning of gene. Please analyse "
                                             "this gene carefully." % gene_id)



                            if summit_peak >= int(feature[1]) and summit_peak <= int(feature[2]): #totally included in feature


                                if feature[0] == "stop_codon" or feature[0] == "start_codon": #if within a one of these codons, associate with exon
                                    logging.info("Summit of the peak %s is located within either start or stop codon. An exon will be associtated to this peak." % fields[0])


                                elif not feature[0] == "tss" and not feature[0] == "tts" and not feature[0] == "CDS": #add feature where summit is included. Exception when the feature is tss or tts. This is associated with the begining of gene
                                    if len(result) == 0:
                                        result.append(feature[0])

                            if feature[0] == "start_codon": #if start codon, check 5' UTR exists and if yes, check if summit is there. If so, replace result with 5' UTR
                                if int(feature[1]) > min_coord and summit_peak < int(feature[1]):
                                    result = []
                                    result.append("5'UTR")

                            elif feature[0] == "stop_codon": #if stop codon, check 3'UTR exists and if yes, check if summit is there.
                                if int(feature[1]) < max_coord and summit_peak > int(feature[2]):
                                    result = []
                                    result.append("3'UTR")

                            other_processed = True



                        #################################################### REVERSE GENES ####################################################
                        if feature[3] == "-":
                            if int(feature[1]) == int(feature[2]) and int(feature[1]) == min_coord and tts == False and other_processed == False:
                                if feature[0] == "tts":
                                    tts = True
                                    ntts +=1
                            elif int(feature[1]) == int(feature[2]) and int(feature[1]) == max_coord and other_processed == True:
                                if feature[0] == "tss":
                                    tss = True
                                    ntss += 1
                            elif int(feature[1]) == int(feature[2]) and int(feature[1]) != min_coord and int(feature[1]) != max_coord :
                                logging.info("Warning. Gene id %s does not have tts or tss starting either at the end or beginning of gene. Please analyse"
                                             "this gene carefully." % gene_id)



                            if summit_peak >= int(feature[1]) and summit_peak <= int(feature[2]): #totally included in feature


                                if feature[0] == "stop_codon" or feature[0] == "start_codon": #if within one of these codons, associate with exon
                                    logging.info("Summit of the peak %s is located within either start or stop codon. An exon will be associtated to this peak." % fields[0])


                                elif not feature[0] == "tss" and not feature[0] == "tts" and not feature[0] == "CDS": #add feature where summit is included. Exception when the feature is tss or tts. This is associated with the begining of gene
                                    if len(result) == 0:
                                        result.append(feature[0])

                            if feature[0] == "start_codon": #if start codon, check 5' UTR exists and if yes, check if summit is there. If so, replace result with 5' UTR
                                if int(feature[1]) < max_coord and summit_peak > int(feature[2]):
                                    result = []
                                    result.append('5\'UTR')

                            elif feature[0] == "stop_codon": #if stop codon, check 3'UTR exists and if yes, check if summit is there.
                                if int(feature[1]) > int(min_coord) and summit_peak < int(feature[1]):
                                    result = []
                                    result.append('3\'UTR')
                            other_processed=True

                    if not result:
                        result =  "['no_feature_predicted']"
                        #logging.info("There is no feature detected in GFF within the gene %s, for which the peak %s was found to be totally within" % (gene_id,fields[0]))

                    fields[i] = fields[i] + ' ' + str(result)
                    outFile.write('\t'.join(fields) + '\n')

            #print(ngenes)
            #print(ntss)
            #print(ntts)




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
