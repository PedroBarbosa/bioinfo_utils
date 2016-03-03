__author__ = 'pedro'

import argparse
import os
import logging
import sys
import csv
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
from collections import defaultdict
import operator
from operator import itemgetter



#def processGTFfile(annotation_dict):

def processFiles(peak_files,threshold, sort,gtf):
    #create file object for the general output file
    if os.path.exists("peaks-general-stats.txt"):
            os.remove("peaks-general-stats.txt")

    #structural annotations variable
    list_ofDicts_annotation = []
    with open(os.getcwd() + "/peaks-general-stats.txt", 'a') as outputFile:

        for filename in peak_files:
            logging.info("Processing %s peaks.." % os.path.basename(filename))

            #general stats
            total_peaks = 0
            peaks_threshold = 0
            scaffolds_peaks = defaultdict(int)
            scaffolds_peaks_threshold = defaultdict(int)
            #output sorted files
            list_of_tuples = []
            #structural annotation
            structural_annotation = defaultdict(tuple)


            #create file object for the sorted peak output file
            if os.path.exists(os.path.splitext(os.path.basename(filename))[0] + "-sortedFoldEnrichment.tsv"):
                os.remove(os.path.splitext(os.path.basename(filename))[0] + "-sortedFoldEnrichment.tsv")
            with open(os.path.splitext(os.path.basename(filename))[0] + "-sortedFoldEnrichment.tsv", "w") as sorted_file:
                writer = csv.writer(sorted_file,dialect=csv.excel_tab)

                with open(filename, 'r') as file:
                    for line_raw in file:
                        line = line_raw.rstrip()
                        if line.strip() and not line.startswith("#") and 'fold_enrichment' in line:
                            header = line.split("\t")
                            if sort:
                                writer.writerow((header))
                        elif line.strip() and not line.startswith("#"):
                            peak = line.split("\t")

                            #general stats
                            total_peaks += 1
                            scaffolds_peaks[peak[0]] += 1

                            #threshold stats
                            if float(peak[7]) >= threshold:
                                peaks_threshold += 1
                                scaffolds_peaks_threshold[peak[0]] += 1

                                #add info about the genome coordinates of the peak
                                peak_without_name = peak[:-1]
                                structural_annotation[peak[-1]] = tuple(peak_without_name)
                                if sort:
                                    list_of_tuples.append((tuple(peak)))

                        elif line.startswith("#") and sort:
                            writer.writerow([line.rstrip()])



 #                   if gtf:
                        #process peaks based on GTF file. Only need the following dict which contains only the peaks above the threshold
#                        processGTFfile(structural_annotation)


                    #append dict for each sample to further process peak annotation
                    #list_ofDicts_annotation.append(structural_annotation)

                    if sort:
                        #sort peaks above threshold by foldEnrichment value and write to the opened file
                        sorted_list = sorted(list_of_tuples,key=itemgetter(7), reverse=True)
                        for processed_peak in sorted_list:
                            writer.writerow((processed_peak))


                    ##############################
                    #Write stats to file##########
                    logging.info("Generating stats..")
                    outputFile.write("########### %s ###########\n" % os.path.basename(filename))
                    outputFile.write("Total number of peaks\t%i\n" % total_peaks)
                    outputFile.write("Number of peaks above the threshold %s\t%i\n\n" % (threshold, peaks_threshold))
                    outputFile.write("Number of different scaffolds/contigs/chromossomes with peaks\t%i\n" % len(scaffolds_peaks))
                    outputFile.write("Average number of peaks per scaffold\t%s\n" % str(round(float(total_peaks)/len(scaffolds_peaks),4)))
                    maxval_1 = max(scaffolds_peaks.iteritems(), key=operator.itemgetter(1))[1]
                    keys_1 = [k for k,v in scaffolds_peaks.items() if v==maxval_1]
                    outputFile.write("Max number of peaks in a scaffold\t%i%s\n\n" % (maxval_1, keys_1))

                    outputFile.write("Number of different scaffolds/contigs/chromossomes with valid peaks above the threshold\t%i\n" % len(scaffolds_peaks_threshold))
                    outputFile.write("Average number of peaks per scaffold considering the peaks above threshold\t%s\n" % str(round(peaks_threshold/len(scaffolds_peaks_threshold),4)))
                    maxval_2 = max(scaffolds_peaks_threshold.iteritems(), key=operator.itemgetter(1))[1]
                    keys_2 = [k for k,v in scaffolds_peaks_threshold.items() if v==maxval_2]
                    outputFile.write("Max number of peaks above the threshold in a scaffold\t%i%s\n\n" % (maxval_2, keys_2))


def main():

    parser = argparse.ArgumentParser(description='Script to analyse fold enrichment in ChipSeq peak files from macs2.')
    parser.add_argument(dest='peak_files', metavar='peakFiles', nargs='+', help='List of files where each line represents a peak file')
    parser.add_argument('--threshold', metavar= 'FLOAT', type=float, required = True, help='Fold enrichment threshold to filter peak files.')
    parser.add_argument('-s', '--sort', action='store_true', help = "Output new filtered files sorted by fold enrichment values above the threshold '-f'. Default: No output files will be written.")
    parser.add_argument('-gff', metavar='--gffFile',nargs=1,help = "Add genome GFF3 file to further analyse the regions of the genome with peaks above the threshold.")
    args = parser.parse_args()


    processFiles(args.peak_files,args.threshold, args.sort, args.gff)


if __name__ == "__main__":
    main()