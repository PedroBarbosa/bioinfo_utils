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
import subprocess
from collections import OrderedDict

def processGFFfile(annotation_dict,gff):

    gff_dict = defaultdict(list)
    no_annotation = 0
    peaks_within_plus = 0
    peaks_within_minus = 0
    partial_peaks_plus = 0
    partial_peaks_minus = 0
    peaks_distance_plus = []
    peaks_distance_minus = []
    peaks2gene_within_plus = {}
    peaks2gene_within_minus = {}
    final_dict = defaultdict(tuple)

    p1 = subprocess.Popen(['cat', gff[0]],stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['grep', '-v', '#'], stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    p3 = subprocess.Popen(['grep','-w','gene'], stdin=p2.stdout, stdout=subprocess.PIPE)
    p2.stdout.close()
    genes = p3.communicate()[0].splitlines()
    p3.stdout.close()

    #add genes to dict
    for gene in genes:
        gene_features = gene.split()
        gff_dict[gene_features[0]].append((gene_features[3], gene_features[4], gene_features[6], gene_features[8]))#scaffold id as key | start, stop, strand and gene_id as list elements as value

    for peak,info in annotation_dict.iteritems():
        scaffold_id = info[0]
        start_peak = info[1]
        length_peak = info[3]
        if scaffold_id in gff_dict:
            genes_per_peak_scaffold = gff_dict[scaffold_id]


            plus_check = True #Flag to avoid searching for genes when the closest by each strand was already found
            minus_check = True
            closer_genes = {}

            for gene in genes_per_peak_scaffold:
                if gene[2] == "+" and plus_check:
                    if start_peak >= gene[0] and start_peak <= gene[1]: #if peak is within a gene in the forward strand
                        peaks_within_plus += 1
                        peaks2gene_within_plus[peak] = gene[3] #dict with the genes that the peaks are within
                        plus_check = False
                        closer_genes['+'] = [gene[3], 'total_within']

                    elif start_peak <= gene[0]: #if peak is upstream a gene in the forward strand, assuming the order of the coordinates in the gff is correct
                        distance_upstream_plus = gene[0] - start_peak
                        if length_peak > distance_upstream_plus:
                            closer_genes['+'] = [gene[3],'partial within '+ [length_peak - distance_upstream_plus]]
                            partial_peaks_plus += 1
                        else:
                            closer_genes['+'] = [gene[3],str(distance_upstream_plus)]
                            peaks_distance_plus.append(distance_upstream_plus)
                        plus_check = False

            for gene in reversed(genes_per_peak_scaffold): #run reversed list to search for minus genes
                if gene[2] == "-" and minus_check:
                    if start_peak >= gene[0] and start_peak <= gene[1]: #if peak is within a gene in the reverse strand
                        peaks_within_minus += 1
                        peaks2gene_within_minus[peak] = gene[3] #dict with the genes that the peaks are within
                        minus_check = False
                        closer_genes['-'] = [gene[3], 'total_within']

                    elif start_peak >= gene[1]:
                        distance_upstream_minus = start_peak - gene[1]
                        if length_peak > distance_upstream_minus:
                            closer_genes['-'] = [gene[3],'partial within '+ [length_peak - distance_upstream_plus]]
                            partial_peaks_minus += 1
                        else:
                            closer_genes['-'] = [gene[3],str(distance_upstream_minus)]
                            peaks_distance_minus.append(distance_upstream_minus)
                        minus_check = False

            if '+' in closer_genes and '-' in closer_genes:
                final_dict[peak] = (info,closer_genes['+'], closer_genes['-'])
            elif '+' in closer_genes:
                final_dict[peak] = (info,closer_genes['+'],'no_gene_minus','-')
            elif '-' in closer_genes:
                final_dict[peak] = (info,'no_gene_plus','-', closer_genes['-'])

        else:
            logging.warning("There are no predicted genes in the scaffold %s where the peak %s was found!" % (scaffold_id,peak))
            no_annotation += 1
            final_dict[peak] = (info,'No annotation available','-','-','-')

    return no_annotation,peaks_within_plus,peaks_within_minus,partial_peaks_plus,partial_peaks_minus,peaks_distance_plus,peaks_distance_minus,peaks2gene_within_plus,\
           peaks2gene_within_minus,final_dict


def processFiles(peak_files,threshold, sort,gff):
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


                    if sort:
                        logging.info("Sorting and writing new peaks file..")
                        #sort peaks above threshold by foldEnrichment value and write to the opened file
                        sorted_list = sorted(list_of_tuples,key=itemgetter(7), reverse=True)
                        for processed_peak in sorted_list:
                            writer.writerow((processed_peak))

                    if gff:
                        logging.info("Processing gff file..")
                        #process peaks based on GTF file. Only need the following dict which contains only the peaks above the threshold
                        no_annotation,peaks_within_plus,peaks_within_minus,partial_peaks_plus,partial_peaks_minus,peaks_distance_plus,peaks_distance_minus,\
                        peaks2gene_within_plus,peaks2gene_within_minus,final_dict = processGFFfile(structural_annotation,gff)


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
                    outputFile.write("Average number of peaks per scaffold considering the peaks above threshold\t%s\n" % str(round(float(peaks_threshold)/len(scaffolds_peaks_threshold),4)))
                    maxval_2 = max(scaffolds_peaks_threshold.iteritems(), key=operator.itemgetter(1))[1]
                    keys_2 = [k for k,v in scaffolds_peaks_threshold.items() if v==maxval_2]
                    outputFile.write("Max number of peaks above the threshold in a scaffold\t%i%s\n\n" % (maxval_2, keys_2))

                    if gff:
                        outputFile.write("##Annotation analysis##\n")
                        outputFile.write("Forward strand:")
                        outputFile.write("\tNumber of peaks totally within the range of genes predicted in the forward strand\t%i\n" % peaks_within_plus)
                        outputFile.write("\tNumber of peaks partially within the range of genes predicted in the forward strand\t%i\n" % partial_peaks_plus)
                        outputFile.write("\tNumber of peaks located upstream of genes predicted in the forward strand\t%i\n" % len(peaks_distance_plus))
                        outputFile.write("\tAverage distance of the upstream peaks to the start of the closest gene predicted in the forward strand\t%s\n" % str(round(float(sum(peaks_distance_plus))/len(peaks_distance_plus),4)))
                        outputFile.write("\tMax distance of a upstream peak to the start of the closest gene predicted in the forward strand\t%i\n" % max(peaks_distance_plus))
                        outputFile.write("\tMin distance of a upstream peak to the start of the closest gene predicted in the forward strand\t%i\n" % min(peaks_distance_plus))
                        outputFile.write("Reverse strand:")
                        outputFile.write("\tNumber of peaks totally within the range of genes predicted in the reverse strand\t%i\n" % peaks_within_minus)
                        outputFile.write("\tNumber of peaks partially within the range of genes predicted in the reverse strand\t%i\n" % partial_peaks_minus)
                        outputFile.write("\tNumber of peaks located upstream of genes predicted in the reverse strand\t%i\n" % len(peaks_distance_minus))
                        outputFile.write("\tAverage distance of the upstream peaks to the start of the closest gene predicted in the reverse strand\t%s\n" % str(round(float(sum(peaks_distance_minus))/len(peaks_distance_minus),4)))
                        outputFile.write("\tMax distance of a upstream peak to the start of the closest gene predicted in the reverse strand\t%i\n" % max(peaks_distance_minus))
                        outputFile.write("\tMin distance of a upstream peak to the start of the closest gene predicted in the reverse strand\t%i\n" % min(peaks_distance_minus))
                        outputFile.write("Number of peaks with no annotation available for the scaffold where they belong\t%i\n\n\n\n" % no_annotation)

                        logging.info("Writing new file with annotation information..")
                        if os.path.exists(os.path.splitext(os.path.basename(filename))[0] + "-annotationInfo.tsv"):
                            os.remove(os.path.splitext(os.path.basename(filename))[0] + "-annotationInfo.tsv")
                        with open(os.path.splitext(os.path.basename(filename))[0] + "-annotationInfo.tsv", "w") as ann_file:
                            for peak, all_info in final_dict:
                                print peak, all_info
            sorted_file.close()
            ann_file.close()
        outputFile.close()


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