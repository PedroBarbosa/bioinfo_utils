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
                    outFile.write("\t".join(line.split("\t")[0:10]) + "\tselected_gene\tfunctional_description\tupstream_distance\tbinding_feature\n")
                elif line.startswith("#") and len(line.split("\t")) == 14:
                    outFile.write("\t".join(line.split("\t")[0:10]) + "\tselected_gene\tupstream_distance\tbinding_feature\n")

                #################################################################################
                #################################################################################
                #################################################################################
                elif "no_gene" in line: #if for one strand there isn't a predicted gene near by, select the other in any case.
                    fields = line.split("\t")
                    if len(fields) == 16: #with annotation
                        upstream_forward = fields[12]
                        upstream_reverse = fields[15]
                        if upstream_reverse == "-" and upstream_forward == "-":
                            logging.error("Error. There is no gene near by the peak %s, when it is supposed to exist at least one. "
                                          "Maybe something was changed in the file manually." % fields[0])
                            exit(1)

                        elif upstream_reverse == "-":
                            if upstream_forward.isdigit():
                                if int(upstream_forward) <= dist:
                                    outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:13]) + '\tpromotor\n')
                                elif int(upstream_forward) > dist:
                                    outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:13]) + '\tintergenic\n')
                            elif "total_within" in upstream_forward:

                                bind_feature = upstream_forward.split("[")[1].replace('\'','').replace('\"','')
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:13]) + '\t' + bind_feature.replace(']','') +'\n')

                            elif "downstream_" in upstream_forward:
                                part_down_distance = int(upstream_forward.split("_")[1])
                                if part_down_distance <= 1000:
                                    outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:13]) + '\t' + 'terminator\n')
                                else:
                                    outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:13]) + '\t' + 'intergenic\n')


                        elif upstream_forward == "-":
                            if upstream_reverse.isdigit():
                                if int(upstream_reverse) <= dist:
                                    outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[13:16]) + '\tpromotor\n')
                                elif int(upstream_reverse) > dist:
                                    outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[13:16]) + '\tintergenic\n')
                            elif "total_within" in upstream_reverse:

                                bind_feature = upstream_reverse.split("[")[1].replace('\'','').replace('\"','')
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[13:16]) + '\t' + bind_feature.replace(']','') +'\n')

                            elif "downstream_" in upstream_reverse:

                                part_down_distance = int(upstream_reverse.split("_")[1])
                                if part_down_distance <= 1000:
                                    outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[13:16]) + '\t' + 'terminator\n')
                                else:
                                    outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[13:16]) + '\t' + 'intergenic\n')

                    elif len(fields) == 14:
                        upstream_forward = fields[11]
                        upstream_reverse = fields[13]
                        if upstream_reverse == "-" and upstream_forward == "-":
                            logging.error("Error. There is no gene near by the peak %s, when it is supposed to exist at least one. "
                                          "Maybe something was changed in the file manually." % fields[0])
                            exit(1)
                    elif upstream_reverse == "-":
                        if upstream_forward.isdigit():
                            if int(upstream_forward) <= dist:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:12]) + '\tpromotor\n')
                            elif int(upstream_forward) > dist:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:12]) + '\tintergenic\n')
                        elif "total_within" in upstream_forward:

                            bind_feature = upstream_forward.split("[")[1].replace('\'','').replace('\"','')
                            outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:12]) + '\t' + bind_feature.replace(']','') +'\n')

                        elif "downstream_" in upstream_forward:
                            part_down_distance = int(upstream_forward.split("_")[1])
                            if part_down_distance <= 1000:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:12]) + '\t' + 'terminator\n')
                            else:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:12]) + '\t' + 'intergenic\n')


                    elif upstream_forward == "-":
                        if upstream_reverse.isdigit():
                            if int(upstream_reverse) <= dist:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[12:14]) + '\tpromotor\n')
                            elif int(upstream_reverse) > dist:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[12:14]) + '\tintergenic\n')
                        elif "total_within" in upstream_reverse:

                            bind_feature = upstream_reverse.split("[")[1].replace('\'','').replace('\"','')
                            outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[12:14]) + '\t' + bind_feature.replace(']','') +'\n')

                        elif "downstream_" in upstream_reverse:

                            part_down_distance = int(upstream_reverse.split(" ")[1])
                            if part_down_distance <= 1000:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[12:14]) + '\t' + 'terminator\n')
                            else:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[12:14]) + '\t' + 'intergenic\n')


                ##########################################################################
                ##########################################################################
                ##########################################################################
                elif not "total_within" in line and not "downstream_" in line \
                     and not "no_gene" in line and not "No annotation available" in line: #if two upstream distances are showed, select the smaller
                    fields = line.split("\t")

                    if len(fields) == 16: #with annotation
                        upstream_forward = int(fields[12])
                        upstream_reverse = int(fields[15])

                        if upstream_forward < upstream_reverse:#select forward gene
                            if upstream_forward <= dist:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:13]) + '\tpromotor\n')
                            else:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:13]) + '\tintergenic\n')


                        elif upstream_forward > upstream_reverse: #select reverse gene

                            if upstream_reverse <= dist:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[13:16]) + '\tpromotor\n')
                            else:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[13:16]) + '\tintergenic\n')

                        else:
                            logging.warning("Peak %s has the same upstream distance for closest forward and reverse gene. Forward gene will be selected.")
                            if upstream_forward <= dist:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[13:16]) + '\tpromotor\n')
                            else:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[13:16]) + '\tintergenic\n')



                    elif len(fields) == 14: #without annotation
                        upstream_forward = int(fields[11])
                        upstream_reverse = int(fields[13])
                        if upstream_forward < upstream_reverse:#select forward gene
                            if upstream_forward <= dist:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:12]) + '\tpromotor\n')
                            else:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:12]) + '\tintergenic\n')


                        elif upstream_forward > upstream_reverse: #select reverse gene
                            if upstream_reverse <= dist:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[12:14]) + '\tpromotor\n')
                            else:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[12:14]) + '\tintergenic\n')

                        else:
                            logging.warning("Peak %s has the same upstream distance for closest forward and reverse gene. Forward gene will be selected." % fields[0])
                            if upstream_forward <= dist:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:12]) + '\tpromotor\n')
                            else:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:12]) + '\tintergenic\n')

                ##########################################################################################
                ##########################################################################################
                ##########################################################################################

                elif "total_within" in line: #if there is one peak totally within a gene, select it.
                    fields = line.split("\t")
                    #indices = [i for i, x in enumerate(fields) if x == "total_within"]

                    indices = [item for item in range(len(fields)) if 'total_within' in fields[item]]
                    if len(indices) > 1:
                        logging.error("Error: Peak %s is totally included in two genes from different strands, this can not happen. Please check if there is any error in the structural "
                                      "annotation file." % fields[0])
                        exit(1)

                    else:
                        i = indices[0]
                        bind_feature = fields[i].split("[")[1].replace('\'','').replace('\"','')
                        if len(fields) == 16:
                            outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[i-2:i+1]) + '\t' + bind_feature.replace(']','') +'\n')
                        elif len(fields) == 14:
                            outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[i-1:i+1]) + '\t' + bind_feature.replace(']','') +'\n')


                #################################################################
                #################################################################
                ################################################################
                elif "downstream_" in line:
                    # 1) If any downstream, select first the upstream gene, if less than dist.
                    # 2) If the upstream is more than dist, select the partial downstream if it exists and has less than 1000bp distance to the end of the gene.
                    # 3) If two partial downstream exist, or one partial and one total, or two total, select the one with less summit distance to the end of the gene
                    # 5) If partial downstream does not exist, select the upstream, even if it's bigger than dist.
                    # 6) If no upstream exists, select the totally downstream gene.


                    fields = line.split("\t")
                    if len(fields) == 16: #with annotation
                        upstream_forward = fields[12]
                        upstream_reverse = fields[15]

                        if upstream_forward.isdigit() and "downstream_" in upstream_reverse and int(upstream_forward) <= dist: #criteria 1
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:13]) + '\tpromotor\n')
                        elif upstream_reverse.isdigit() and "downstream_" in upstream_forward and int(upstream_reverse) <= dist: #criteria 1
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[13:16]) + '\tpromotor\n')

                        elif upstream_reverse.isdigit() and int(upstream_reverse) > dist and "downstream_" in upstream_forward: #criteria 2
                            if int(upstream_forward.split("_")[1]) >= 1000:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[13:16]) + '\tintergenic\n')
                            else:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:13]) + '\tterminator\n')

                        elif upstream_forward.isdigit() and int(upstream_forward) > dist and "downstream_" in upstream_reverse: #criteria 2
                            if int(upstream_reverse.split("_")[1]) >= 1000:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:13]) + '\tintergenic\n')
                            else:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[13:16]) + '\tterminator\n')

                        elif "downstream_" in upstream_reverse and "downstream_" in upstream_forward: #criteria 3
                            dist_down_for = int(upstream_forward.split("_")[1])
                            dist_down_rev = int(upstream_reverse.split("_")[1])
                            if dist_down_for <= dist_down_rev:
                                if dist_down_for <= 1000:
                                    outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:13]) + '\tterminator\n')
                                else:
                                    outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:13]) + '\tintergenic\n')

                            else:
                                if dist_down_rev <= 1000:
                                    outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[13:16]) + '\tterminator\n')
                                else:
                                    outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[13:16]) + '\tintergenic\n')


   #                     elif "total_" in upstream_forward: #criteria 6
   #                         outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:13]) + '\tintergenic\n')

   #                     elif "total_" in upstream_reverse:
   #                         outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[13:16]) + '\tintergenic\n')


                    elif len(fields) == 14: #without annotation

                        upstream_forward = fields[11]
                        upstream_reverse = fields[13]

                        if upstream_forward.isdigit() and "downstream_" in upstream_reverse and int(upstream_forward) <= dist: #criteria 1
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:12]) + '\tpromotor\n')
                        elif upstream_reverse.isdigit() and "downstream_" in upstream_forward and int(upstream_reverse) <= dist: #criteria 1
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[12:14]) + '\tpromotor\n')

                        elif upstream_reverse.isdigit() and int(upstream_reverse) > dist and "downstream_" in upstream_forward: #criteria 2
                            if int(upstream_forward.split("_")[1]) >= 1000:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[12:14]) + '\tintergenic\n')
                            else:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:12]) + '\tterminator\n')

                        elif upstream_forward.isdigit() and int(upstream_forward) > dist and "downstream_" in upstream_reverse: #criteria 2
                            if int(upstream_reverse.split("_")[1]) >= 1000:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:12]) + '\tintergenic\n')
                            else:
                                outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[12:14]) + '\tterminator\n')

                        elif "downstream_" in upstream_reverse and "downstream_" in upstream_forward: #criteria 3
                            dist_down_for = int(upstream_forward.split("_")[1])
                            dist_down_rev = int(upstream_reverse.split("_")[1])
                            if dist_down_for <= dist_down_rev:
                                if dist_down_for <= 1000:
                                    outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:12]) + '\tterminator\n')
                                else:
                                    outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:12]) + '\tintergenic\n')

                            else:
                                if dist_down_rev <= 1000:
                                    outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[12:14]) + '\tterminator\n')
                                else:
                                    outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[12:14]) + '\tintergenic\n')


        #                elif "total_" in upstream_forward: #criteria 6
        #                    outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[10:12]) + '\tintergenic\n')

        #                elif "total_" in upstream_reverse:
        #                    outFile.write("\t".join(fields[0:10]) + "\t" + '\t'.join(fields[12:14]) + '\tintergenic\n')

                else:
                    fields = line.split("\t")
                    logging.warning("Peak %s has no gene associated in any of the strands. This scaffold does not contain any predicted gene in "
                                    "the annotation file." % fields[0])
                    #logging.warning(line)
                    if len(fields) == 16:
                        outFile.write("\t".join(fields[0:10]) + "\tNo annotation available\t-\t-\tunknown\n")
                    elif len(fields) == 14:
                        outFile.write("\t".join(fields[0:10]) + "\tNo annotation available\t-\tunknown\n")





def main():
    parser = argparse.ArgumentParser(description='Script to choose the closest gene of each peak (reverse or forward) based on some annotation criteria. Requires the output of '
                'getPeakLocationWithinGeneBoundariesFromSummit script as the input file. Criteria to select the peaks:'
                '1- If for one strand there isn\'t a predicted gene near by, select the other in any case.'
                '2 - If the peak summit is totally within gene boundaries (gene boundarie is defined as the begining of the genes, which '
                'does not imply to be the transcription initiation site (tss).'
                '3 - If the peak summit is upstream the beginning of gene less than "--upstream" distance.'
                '4 - If the peak summit is downstream of the end of gene less than 1000bp.')


    parser.add_argument(dest='peak_file', metavar='peakFile', nargs=1,
                        help='Excel file with information about peak location regarding the closest gene forward and reverse.')
    parser.add_argument("-u", metavar="--upstream", type=int, required=True,help="Maximum uptream distance of the peak to a gene to be considered promotor region." )
    parser.add_argument("-o", metavar='--outputFile', required=True, help="Output file to write the new excel.")
    args = parser.parse_args()

    processPeaksFile(args.peak_file,args.u, args.o)

if __name__ == "__main__":
    main()