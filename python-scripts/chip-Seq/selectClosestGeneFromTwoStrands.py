import argparse

def processPeaksFile(peakFile,dist, outputFile):

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

                            if int(feature[1]) >= int(start_peak) and int(feature[1]) <= int(end_peak):
                                result.append(feature[0])

                        if int(start_peak) >= int(feature[1]) and int(start_peak) <= int(feature[2]): #totally included in feature
                            result.append(feature[0])


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
                        result =  '[No feature predicted]'
                        logging.info("There is no feature detected in GFF within the gene %s, for which the peak %s was found to be totally within" % (gene_id,fields[0]))
                    fields[i] = fields[i] + ' ' + str(result)
                    outFile.write('\t'.join(fields) + '\n')
                    #print(fields)

def main():
    parser = argparse.ArgumentParser(description='Script to choose the closest gene of each peak (reverse or forward) based on some annotation criteria. Requires the output of '
                                                 'getPeakLocationWithinGeneBoundaries script as the input file. Criteria to select the peaks:'
                                                 '1 - If the gene is totally within gene boundaries (.'
                                                 '2 - If the gene is partially wihtin gene boundaries')
    parser.add_argument(dest='peak_file', metavar='peakFile', nargs=1,
                        help='Excel file with information about peak location regarding the closest gene forward and reverse.')
    parser.add_argument("-u", metavar="--upstream", type=int, required=True, nargs=1,help="Maximum uptream distance of the peak to a gene to be considered promotor region." )
    parser.add_argument("-o", metavar='--outputFile', required=True, help="Output file to write the new excel.")
    args = parser.parse_args()

    processPeaksFile(args.peak_file,args.u, args.o)

if __name__ == "__main__":
    main()