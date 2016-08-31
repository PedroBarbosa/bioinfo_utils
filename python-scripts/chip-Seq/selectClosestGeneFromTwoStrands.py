import argparse
import logging
import sys
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')

def processPeaksFile(peakFile,dist, outputFile):

    with open(peakFile[0], 'r') as peakFile:
        with open(outputFile, 'w') as outFile:
            for line in peakFile:
                line =line.rstrip()
                if line.startswith("#") and len(line.split("\t")) == 16:
                    outFile.write("\t".join(line.split("\t")[0:9]) + "\tselected_gene\tfunctional_description\tupstream_distance\tbinding_feature\n")
                elif line.startswith("#") and len(line.split("\t")) == 14:
                    outFile.write("\t".join(line.split("\t")[0:9]) + "\tselected_gene\tupstream_distance\tbinding_feature\n")

                #################################################################################
                #################################################################################
                #################################################################################
                elif not "total" in line and not "partial" and not "no_gene":#if two upstream distances are showed, select the smaller
                    fields = line.split("\t")

                    if len(fields) == 16: #with annotation
                        upstream_forward = int(fields[12])
                        upstream_reverse = int(fields[15])
                        if upstream_forward < upstream_reverse:#select forward gene
                            if upstream_forward <= dist:
                                outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[10:12]) + '\tpromotor')
                            else:
                                outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[10:12]) + '\tintergenic')


                        elif upstream_forward > upstream_reverse: #select reverse gene
                            if upstream_reverse <= dist:
                                outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[13:15]) + '\tpromotor')
                            else:
                                outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[13:15]) + '\tintergenic')

                        else:
                            logging.warning("Peak %s has the same upstream distance for closest forward and reverse gene. Forward gene will be selected.")
                            if upstream_forward <= dist:
                                outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[13:15]) + '\tpromotor')
                            else:
                                outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[13:15]) + '\tintergenic')



                    elif len(fields) == 14: #without annotation
                        upstream_forward = int(fields[11])
                        upstream_reverse = int(fields[13])
                        if upstream_forward < upstream_reverse:#select forward gene
                            if upstream_forward <= dist:
                                outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[10:11]) + '\tpromotor')
                            else:
                                outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[10:11]) + '\tintergenic')


                        elif upstream_forward > upstream_reverse: #select reverse gene
                            if upstream_reverse <= dist:
                                outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[12:13]) + '\tpromotor')
                            else:
                                outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[12:13]) + '\tintergenic')

                        else:
                            logging.warning("Peak %s has the same upstream distance for closest forward and reverse gene. Forward gene will be selected." % fields[0])
                            if upstream_forward <= dist:
                                outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[10:11]) + '\tpromotor')
                            else:
                                outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[10:11]) + '\tintergenic')

                ##########################################################################################
                ##########################################################################################
                ##########################################################################################
                elif "total_within" in line: #if there is on peak totally within a gene, select it.
                    fields = line.split("\t")
                    indices = [i for i, x in enumerate(fields) if x == "total_within"]
                    if len(indices) > 1:
                        logging.error("Error: Peak %s is totally included in two genes from different strands, this can not happen. Please check if there is any error in the structural "
                                      "annotation file." % fields[0])
                        exit(1)
                    elif "partial" in line:
                        logging.error("Error: Peak %s is tottally included in one gene and partialle included in other in the opposite strand. Please check for errors"
                                      "in the annotation file." % fields[0])
                        exit(1)

                    else:
                        i = indices[0]
                        if len(fields) == 16:
                            bind_feature = fields[i].split("[")[1].replace('\'','')
                            outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[i-2:i]) + '\t' + bind_feature.replace(']',''))
                        elif len(fields) == 14:
                            bind_feature = fields[i].split("[")[1].replace('\'','')
                            outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[i-1:i]) + '\t' + bind_feature.replace(']',''))


                ##############################################################
                ##############################################################
                ##############################################################
                elif "partial within" in line: # if there isn't any totally within, select the upstream partial ones.
                    fields = line.split("\t")
                    indices = [i for i, x in enumerate(fields) if "partial within" in x]
                    if len(indices) > 1:
                        logging.error("Error: Peak %s is partially included in two genes from different strands, this can not happen. Please check if there is any error in the structural "
                                      "annotation file." % fields[0])
                        exit(1)
                    else:
                        i = indices[0]
                        if len(fields) == 16:
                            outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[i-2:i]) + '\tpartial_upstream')
                        elif len(fields) == 14:
                            outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[i-1:i]) + '\tpartial_upstream')



                #################################################################
                #################################################################
                ################################################################
                elif "downstream" in line:
                    # 1) If any downstream, select first the upstream gene, if less than dist.
                    # 2) If the upstream is more than dist, select the partial downstream if it exists.
                    # 3) If partial downstream does not exist, select the upstream, even if it's bigger than dist.
                    # 4) If no upstream exists, select the totally downstream gene.

                    if len(fields) == 16: #with annotation
                        upstream_forward = fields[12]
                        upstream_reverse = fields[15]

                        if upstream_forward.isdigit():
                            if not "downstream" in upstream_forward and int(upstream_forward) <= dist: #criteria 1
                                outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[10:12]) + '\tpromotor')
                        elif upstream_reverse.isdigit():
                            if not "downstream" in upstream_reverse and int(upstream_reverse) <= dist: #criteria 1
                                outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[13:15]) + '\tpromotor')

                        if upstream_reverse.isdigit():
                            if "partial" in upstream_forward and int(upstream_reverse) > dist: #criteria 2
                                outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[10:12]) + '\tpartial_downstream')

                        elif upstream_forward.isdigit():
                            if "partial" in upstream_reverse and int(upstream_forward) > dist: #criteria 2
                                outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[13:15]) + '\tpartial_downstream')

                        if upstream_forward.isdigit() and "totally" in upstream_reverse:
                            outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[10:12]) + '\tintergenic')

                        elif upstream_reverse.isdigit() and "totally" in upstream_forward: #criteria 3
                            outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[13:15]) + '\tintergenic') #criteria 3

                        elif "totally" in upstream_forward and "totally" in upstream_reverse:
                            logging.warning("Something wrong is happening. Peak %s can't be located downstream of two genes." % fields[0])

                        elif "totally" in upstream_forward:
                            outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[10:12]) + '\tintergenic')

                        elif "totally" in upstream_reverse:
                            outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[13:15]) + '\tintergenic')


                    elif len(fields) == 14: #without annotation

                        upstream_forward = fields[11]
                        upstream_reverse = fields[13]

                        if upstream_forward.isdigit():
                            if not "downstream" in upstream_forward and int(upstream_forward) <= dist: #criteria 1
                                outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[10:11]) + '\tpromotor')
                        elif upstream_reverse.isdigit():
                            if not "downstream" in upstream_reverse and int(upstream_reverse) <= dist: #criteria 1
                                outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[12:13]) + '\tpromotor')

                        if upstream_reverse.isdigit():
                            if "partial" in upstream_forward and int(upstream_reverse) > dist: #criteria 2
                                outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[10:11]) + '\tpartial_downstream')

                        elif upstream_forward.isdigit():
                            if "partial" in upstream_reverse and int(upstream_forward) > dist: #criteria 2
                                outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[12:13]) + '\tpartial_downstream')

                        if upstream_forward.isdigit() and "totally" in upstream_reverse:
                            outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[10:11]) + '\tintergenic')

                        elif upstream_reverse.isdigit() and "totally" in upstream_forward: #criteria 3
                            outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[12:13]) + '\tintergenic') #criteria 3

                        elif "totally" in upstream_forward and "totally" in upstream_reverse:
                            logging.warning("Something wrong is happening. Peak %s can't be located downstream of two genes." % fields[0])

                        elif "totally" in upstream_forward:
                            outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[10:11]) + '\tintergenic')

                        elif "totally" in upstream_reverse:
                            outFile.write("\t".join(fields[0:9]) + "\t" + '\t'.join(fields[12:13]) + '\tintergenic')


                else:
                    logging.warning("Peak %s has no gene associated in any of the strands. It's likely that this scaffold does not contain any predicted gene in "
                                    "the annotation file." % fields[0])
                    logging.warning(line)




def main():
    parser = argparse.ArgumentParser(description='Script to choose the closest gene of each peak (reverse or forward) based on some annotation criteria. Requires the output of '
                                                 'getPeakLocationWithinGeneBoundaries script as the input file. Criteria to select the peaks:'
                                                 '1 - If the gene is totally within gene boundaries (gene boundarie is defined as the begining onf the genes, which '
                                                 'does not imply to always be the transcription initiation site (tss).'
                                                 '2 - If the gene is partially within gene boundaries.'
                                                 '3 - If the peak is upstream the gene less than "--upstream" distance.')

    parser.add_argument(dest='peak_file', metavar='peakFile', nargs=1,
                        help='Excel file with information about peak location regarding the closest gene forward and reverse.')
    parser.add_argument("-u", metavar="--upstream", type=int, required=True, nargs=1,help="Maximum uptream distance of the peak to a gene to be considered promotor region." )
    parser.add_argument("-o", metavar='--outputFile', required=True, help="Output file to write the new excel.")
    args = parser.parse_args()

    processPeaksFile(args.peak_file,args.u, args.o)

if __name__ == "__main__":
    main()