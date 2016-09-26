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


def processGFFfile(annotation_dict,gff):

    gff_dict = defaultdict(list)
    no_annotation = 0
    peaks_downstream_plus = 0
    peaks_downstream_minus = 0
    peaks_partially_downstream_plus = 0
    peaks_partially_downstream_minus = 0
    peaks_within_plus = 0
    peaks_within_minus = 0
    partial_peaks_plus = 0
    partial_peaks_minus = 0
    peaks_distance_plus = []
    peaks_distance_minus = []
    peaks2gene_within_plus = {}
    peaks2gene_within_minus = {}
    final_dict = defaultdict(tuple)
    no_gene_plus = 0
    no_gene_minus = 0
    count=0
    peak_count = 0


    p1 = subprocess.Popen(['cat', gff[0]],stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['grep', '-v', '#'], stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    p3 = subprocess.Popen(['grep','-w','gene'], stdin=p2.stdout, stdout=subprocess.PIPE, universal_newlines=True)
    p2.stdout.close()
    genes = p3.communicate()[0].splitlines()
    p3.stdout.close()

    logging.info("Gff read.")
    #add genes to dict
    for gene in genes:

        gene_features = gene.split()
        gff_dict[gene_features[0]].append((gene_features[3], gene_features[4], gene_features[6], gene_features[8]))#scaffold id as key | start, stop, strand and gene_id as list elements as value

  #  for k in iter(gff_dict.keys()):
  #      print(k)
    for peak,info in iter(annotation_dict.items()):

        peak_count +=1

        scaffold_id = info[0]
        peak_summit = int(info[4])
        if scaffold_id in gff_dict:

            genes_per_peak_scaffold = gff_dict[scaffold_id]
            plus_check = True #Flag to avoid searching for genes when the closest by each strand was already found
            minus_check = True
            exist_gene_plus = False
            exist_gene_minus = False
            temp_downstream_dist_plus = 0
            temp_downstream_dist_minus = 0
            closer_genes = {}
            temp_plus = []
            temp_minus = []


            for gene in genes_per_peak_scaffold:
                start = int(gene[0])
                end = int(gene[1])
                if gene[2] == "+" and plus_check:

                    exist_gene_plus = True
                    if peak_summit >= start and peak_summit <= end:
                        #if start_peak >= start and start_peak <= end and end_peak <= end: #if peak is within a gene in the forward strand
                        peaks_within_plus += 1
                        peaks2gene_within_plus[peak] = gene[3] #dict with the genes that the peaks are within
                        plus_check = False
                        closer_genes['+'] = [gene[3], 'total_within']
                        break
                    #elif start_peak <= start: #if peak is upstream a gene in the forward strand, assuming the order of the coordinates in the gff is correct
                    elif peak_summit < start:

                        distance_upstream_plus = start - peak_summit
                        if len(temp_plus) != 0:
                            if temp_downstream_dist_plus <= 1000 and distance_upstream_plus > 2000:
                                #print(temp_plus)
                                closer_genes['+'] = temp_plus
                                peaks_partially_downstream_plus += 1
                                plus_check = False
                                break

                            else:
                                closer_genes['+'] = [gene[3],str(distance_upstream_plus)]
                                peaks_distance_plus.append(distance_upstream_plus)
                                plus_check = False
                                break
                        else:

                            closer_genes['+'] = [gene[3],str(distance_upstream_plus)]
                            peaks_distance_plus.append(distance_upstream_plus)
                            plus_check = False
                            break

                    elif peak_summit > end and plus_check:#if peak is partially downtream all genes in scaffold
                        #temp_plus = [gene[3], 'partial_downstream ' +str(end - start_peak)]

                        temp_plus = [gene[3], 'downstream_' + str(peak_summit - end)]
                        temp_downstream_dist_plus = peak_summit - end



            if plus_check and exist_gene_plus: #checked all genes for upstream location in the forward strand but didn't found any, despite the fact it found downstream

                peaks_downstream_plus += 1
                closer_genes['+'] = temp_plus

            elif exist_gene_plus == False:
                no_gene_plus += 1


            #######################################################
            ######################################################
            #####################################################
            for gene in reversed(genes_per_peak_scaffold): #run reversed list to search for minus genes

                end = int(gene[0])
                start = int(gene[1])
                if gene[2] == "-" and minus_check:

                    exist_gene_minus = True
                    if peak_summit <= start and peak_summit >= end:
                        #if start_peak >= start and start_peak <= end and end_peak <= end: #if peak is within a gene in the forward strand
                        peaks_within_minus += 1
                        peaks2gene_within_minus[peak] = gene[3] #dict with the genes that the peaks are within
                        minus_check = False
                        closer_genes['-'] = [gene[3], 'total_within']
                        break
                    #elif start_peak <= start: #if peak is upstream a gene in the forward strand, assuming the order of the coordinates in the gff is correct
                    elif peak_summit > start:

                        distance_upstream_minus = peak_summit - start
                        if len(temp_minus) != 0:
                            if temp_downstream_dist_minus <= 1000 and distance_upstream_minus > 2000:
                                #print(temp_plus)
                                closer_genes['-'] = temp_minus
                                peaks_partially_downstream_minus += 1
                                minus_check = False
                                break

                            else:
                                closer_genes['-'] = [gene[3],str(distance_upstream_minus)]
                                peaks_distance_minus.append(distance_upstream_minus)
                                minus_check = False
                                break
                        else:

                            closer_genes['-'] = [gene[3],str(distance_upstream_minus)]
                            peaks_distance_minus.append(distance_upstream_minus)
                            minus_check = False
                            break

                    elif peak_summit < end and minus_check:#if peak is partially downtream all genes in scaffold

                        temp_minus = [gene[3], 'downstream_' + str(end - peak_summit)]
                        temp_downstream_dist_minus = end - peak_summit



            if minus_check and exist_gene_minus: #checked all genes for upstream location in the forward strand but didn't found any, despite the fact it found downstream

                peaks_downstream_minus += 1
                closer_genes['-'] = temp_minus

            elif exist_gene_minus == False:
                no_gene_minus += 1



            if '+' in closer_genes and '-' in closer_genes:
                final_dict[peak] = info + (closer_genes['+'][0],closer_genes['+'][1], closer_genes['-'][0], closer_genes['-'][1],)
            elif '+' in closer_genes:
                #print(info)
                #print(closer_genes['+'])
                final_dict[peak] = info + (closer_genes['+'][0], closer_genes['+'][1],'no_gene_minus','-',)
            elif '-' in closer_genes:
                final_dict[peak] = info + ('no_gene_plus','-', closer_genes['-'][0], closer_genes['-'][1],)



        else:
            #logging.warning("There are no predicted genes in the scaffold %s where the peak %s was found!" % (scaffold_id,peak))
            no_annotation += 1
            final_dict[peak] = info + ('No annotation available','-','-','-',)



    #print(no_annotation,peaks_within_plus,peaks_within_minus,partial_peaks_plus,partial_peaks_minus,peaks_distance_plus,peaks_distance_minus,peaks2gene_within_plus,\
     #      peaks2gene_within_minus,peaks_downstream_plus, peaks_partially_downstream_plus, peaks_downstream_minus, peaks_partially_downstream_minus, final_dict)
    return no_annotation,peaks_within_plus,peaks_within_minus,partial_peaks_plus,partial_peaks_minus,peaks_distance_plus,peaks_distance_minus,peaks2gene_within_plus,\
           peaks2gene_within_minus,peaks_downstream_plus, peaks_partially_downstream_plus, peaks_downstream_minus, peaks_partially_downstream_minus, final_dict, no_gene_plus, no_gene_minus




def processFiles(peak_files,threshold, sort,gff, addAnnotation, col):
    #create file object for the general output file
    if os.path.exists("peaks-general-stats.txt"):
            os.remove("peaks-general-stats.txt")

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
                        sorted_list = sorted(list_of_tuples, key=lambda x: float(x[7]), reverse=True)
                        #sorted_list = sorted(list_of_tuples,key=itemgetter(7), reverse=True) doesnt work for doubles bigger than 10
                        for processed_peak in sorted_list:
                            writer.writerow((processed_peak))
                    else:
                        os.remove(os.path.splitext(os.path.basename(filename))[0] + "-sortedFoldEnrichment.tsv")
                    if gff:
                        logging.info("Processing gff file..")
                        #process peaks based on GTF file. Only need the following dict which contains only the peaks above the threshold
                        no_annotation,peaks_within_plus,peaks_within_minus,partial_peaks_plus,partial_peaks_minus,peaks_distance_plus,peaks_distance_minus,\
                        peaks2gene_within_plus,peaks2gene_within_minus,peaks_downstream_plus, peaks_partially_downstream_plus, peaks_downstream_minus,\
                        peaks_partially_downstream_minus,final_dict,no_gene_plus, no_gene_minus = processGFFfile(structural_annotation,gff)



                    ##############################
                    #Write stats to file##########

                    logging.info("Generating stats..")
                    outputFile.write("########### %s ###########\n" % os.path.basename(filename))

                    if total_peaks == 0:
                        logging.info("No peaks detected in file|")
                        outputFile.write("No peaks in file!!\n\n\n\n\n\n")

                    else:
                        outputFile.write("Total number of peaks\t%i\n" % total_peaks)
                        outputFile.write("Number of peaks above the threshold %s\t%i\n\n" % (threshold, peaks_threshold))
                        outputFile.write("Number of different scaffolds/contigs/chromossomes with peaks\t%i\n" % len(scaffolds_peaks))
                        outputFile.write("Average number of peaks per scaffold\t%s\n" % str(round(float(total_peaks)/len(scaffolds_peaks),4)))
                        maxval_1 = max(iter(scaffolds_peaks.items()), key=operator.itemgetter(1))[1]
                        keys_1 = [k for k,v in iter(scaffolds_peaks.items()) if v==maxval_1]
                        outputFile.write("Max number of peaks in a scaffold\t%i%s\n\n" % (maxval_1, keys_1))

                        if len(scaffolds_peaks_threshold) > 0:
                            outputFile.write("Number of different scaffolds/contigs/chromossomes with valid peaks above the threshold\t%i\n" % len(scaffolds_peaks_threshold))
                            outputFile.write("Average number of peaks per scaffold considering the peaks above threshold\t%s\n" % str(round(float(peaks_threshold)/len(scaffolds_peaks_threshold),4)))
                            maxval_2 = max(iter(scaffolds_peaks_threshold.items()), key=operator.itemgetter(1))[1]
                            keys_2 = [k for k,v in iter(scaffolds_peaks_threshold.items()) if v==maxval_2]
                            outputFile.write("Max number of peaks above the threshold in a scaffold\t%i%s\n\n" % (maxval_2, keys_2))


                            if gff:

                                outputFile.write("##Annotation analysis of the peaks above the threshold:##\n")
                                outputFile.write("Forward strand:\n")
                                outputFile.write("\tNumber of peaks totally within the range of genes predicted in the forward strand\t%i\n" % peaks_within_plus)
                                #outputFile.write("\tNumber of peaks partially within the range of genes predicted in the forward strand\t%i\n" % partial_peaks_plus)
                                outputFile.write("\tNumber of peaks located upstream of genes predicted in the forward strand\t%i\n" % len(peaks_distance_plus))
                                outputFile.write("\tNumber of peaks where all the genes predicted in the forward strand are totally upstream of the peak itself (peaks located downstream of all genes of the scaffold)\t%i\n" % peaks_downstream_plus)
                                outputFile.write("\tNumber of peaks where located downstream of one gene up to 1000bp, considering that there is a gene more than 2000bp downstream of the peak)\t%i\n" % peaks_partially_downstream_plus)
                                #outputFile.write("\tNumber of peaks with no gene information in the forward strand (scaffold of the peak has only gene predictions in reverse strand)\t%i\n" % int(peaks_threshold - sum([peaks_within_plus,len(peaks_distance_plus),peaks_downstream_plus,peaks_partially_downstream_plus,no_annotation])))
                                outputFile.write("\tNumber of peaks with no gene information in the forward strand (scaffold of the peak has only gene predictions in reverse strand)\t%i\n" % no_gene_plus)

                                if len(peaks_distance_plus) > 0:
                                    outputFile.write("\tAverage distance of the upstream peaks to the start of the closest gene predicted in the forward strand\t%s\n" % str(round(float(sum(peaks_distance_plus))/len(peaks_distance_plus),4)))
                                    outputFile.write("\tMax distance of an upstream peak to the start of the closest gene predicted in the forward strand\t%i\n" % max(peaks_distance_plus))
                                    outputFile.write("\tMin distance of an upstream peak to the start of the closest gene predicted in the forward strand\t%i\n" % min(peaks_distance_plus))
                                else:
                                    outputFile.write("\tAverage distance of the upstream peaks to the start of the closest gene predicted in the forward strand\t%s\n" % str(0))
                                    outputFile.write("\tMax distance of an upstream peak to the start of the closest gene predicted in the forward strand\t%i\n" % 0)
                                    outputFile.write("\tMin distance of an upstream peak to the start of the closest gene predicted in the forward strand\t%i\n" % 0)

                                outputFile.write("Reverse strand:\n")
                                outputFile.write("\tNumber of peaks totally within the range of genes predicted in the reverse strand\t%i\n" % peaks_within_minus)
                                #outputFile.write("\tNumber of peaks partially within the range of genes predicted in the reverse strand\t%i\n" % partial_peaks_minus)
                                outputFile.write("\tNumber of peaks located upstream of genes predicted in the reverse strand\t%i\n" % len(peaks_distance_minus))
                                outputFile.write("\tNumber of peaks where all the genes predicted in the reverse strand are totally upstream of the peak itself (peaks located downstream of all genes of the scaffold)\t%i\n" % peaks_downstream_minus)
                                outputFile.write("\tNumber of peaks where located downstream of one gene up to 1000bp, considering that there is a gene more than 2000bp downstream of the peak)\t%i\n" % peaks_partially_downstream_minus)
                                #outputFile.write("\tNumber of peaks with no gene information in the reverse strand (scaffold of the peak has only gene predictions in forward strand)\t%i\n" % int(peaks_threshold - sum([peaks_within_minus,len(peaks_distance_minus),peaks_downstream_minus,peaks_partially_downstream_minus, no_annotation])))
                                outputFile.write("\tNumber of peaks with no gene information in the reverse strand (scaffold of the peak has only gene predictions in forward strand)\t%i\n" % no_gene_minus)
                                if len(peaks_distance_minus) > 0:
                                    outputFile.write("\tAverage distance of the upstream peaks to the start of the closest gene predicted in the reverse strand\t%s\n" % str(round(float(sum(peaks_distance_minus))/len(peaks_distance_minus),4)))
                                    outputFile.write("\tMax distance of an upstream peak to the start of the closest gene predicted in the reverse strand\t%i\n" % max(peaks_distance_minus))
                                    outputFile.write("\tMin distance of an upstream peak to the start of the closest gene predicted in the reverse strand\t%i\n" % min(peaks_distance_minus))
                                else:
                                    outputFile.write("\tAverage distance of the upstream peaks to the start of the closest gene predicted in the reverse strand\t%s\n" % str(0))
                                    outputFile.write("\tMax distance of an upstream peak to the start of the closest gene predicted in the reverse strand\t%i\n" % 0)
                                    outputFile.write("\tMin distance of an upstream peak to the start of the closest gene predicted in the reverse strand\t%i\n" % 0)
                                outputFile.write("Number of peaks with no annotation available for the scaffold where they belong\t%i\n\n\n\n\n\n" % no_annotation)

                                logging.info("Writing new file with annotation information..")
                                if os.path.exists(os.path.splitext(os.path.basename(filename))[0] + "-annotationInfo.tsv"):
                                    os.remove(os.path.splitext(os.path.basename(filename))[0] + "-annotationInfo.tsv")
                                with open(os.path.splitext(os.path.basename(filename))[0] + "-annotationInfo.tsv", "w") as ann_file:
                                    writer_ann = csv.writer(ann_file,dialect=csv.excel_tab)
                                    if addAnnotation:
                                        writer_ann.writerow(["#peak_name","#scaffold_id","#start","#end","#length","#abs_summit","#pileup","#-log10(pvalue)","#fold_enrichment", \
                                    "#-log10(qvalue)", "#closest_gene_forward", "#functional_description_forward","#upstream_dist_forward","#closest_gene_reverse","#functional_description_reverse"\
                                    ,"#upstream_dist_reverse"])

                                        final_dict_updated = processFromBlastTab(final_dict,addAnnotation,col)
                                        for peak, all_info in iter(final_dict_updated.items()):
                                            #info = "\t".join(all_info).replace("\"","")
                                            writer_ann.writerow((peak,'\t'.join(all_info)))
                                    else:
                                        #sorted_dict = sorted(final_dict..items(),key=lambda (k,v): v(8),reverse=True)
                                        writer_ann.writerow(["#peak_name","#scaffold_id","#start","#end","#length","#abs_summit","#pileup","#-log10(pvalue)","#fold_enrichment", \
                                                             "#-log10(qvalue)", "#closest_gene_forward", "#upstream_dist_forward","#closest_gene_reverse","#upstream_dist_reverse"])
                                        for peak, all_info in iter(final_dict.items()):
                                            #info = "\t".join(all_info).replace("\"","")
                                            writer_ann.writerow((peak,'\t'.join(all_info)))
                                ann_file.close()
                                removeChar(filename)
                        else:
                            logging.info("No peaks detected above the threshold.")
                            outputFile.write("No peaks detected above the threshold!\n\n\n\n\n\n")


        sorted_file.close()
        outputFile.close()


def processFromBlastTab(final_dict, functionalAnnotation,col):
    dict_funct={}
    annotated_features = 0
    logging.info("Processing annotation file ..")
    with open(functionalAnnotation) as file:
        previous_query = ""

        for line in file:
            line.rstrip()
            if not line.startswith('#') :
                query = line.split()[0]
                hit = ""
                if '.' in query:
                    query = query.split('.')[0]

                if not query in dict_funct:

                    if query == previous_query:
                        previous_query = query
                    else:
                        annotated_features += 1
                        for c in col:
                            hit += line.split('\t')[c-1] + " "
                        dict_funct[query] = hit


    for peak,all_info in iter(final_dict.items()):
        list_info = list(all_info)
        gene_forward = list_info[9]
        gene_reverse = list_info[11]

        #no structural annotation, thus no functional annotation
        if gene_forward == "No annotation available":
                list_info.insert(10,'-')
                list_info.insert(13,'-')
                final_dict[peak] = tuple(list_info)


        #peaks where genes in both strands were identified
        elif gene_forward != "no_gene_plus" and gene_reverse != "no_gene_minus":
            #there is functional annotation for the genes predicted in plus and minus strand
            if gene_forward in dict_funct and gene_reverse in dict_funct:
                list_info.insert(10,dict_funct[gene_forward])
                list_info.insert(13,dict_funct[gene_reverse])
                final_dict[peak] = tuple(list_info)

            elif gene_forward in dict_funct:
                list_info.insert(10,dict_funct[gene_forward])
                list_info.insert(13,'No function detected')
                final_dict[peak] = tuple(list_info)
            elif gene_reverse in dict_funct:
                list_info.insert(10,'No function detected')
                list_info.insert(13,dict_funct[gene_reverse])
                final_dict[peak] = tuple(list_info)
            else:
                list_info.insert(10,'No function detected')
                list_info.insert(13,'No function detected')
                final_dict[peak] = tuple(list_info)

        #peaks with only genes in reverse strand
        elif gene_forward == "no_gene_plus" and gene_reverse != "no_gene_minus":
            #peaks where the reverse gene has functional annotation
            if gene_reverse in dict_funct:
                list_info.insert(10,'-')
                list_info.insert(13,dict_funct[gene_reverse])
                final_dict[peak] = tuple(list_info)
            else:
                list_info.insert(10,'-')
                list_info.insert(13,'No function detected')
                final_dict[peak] = tuple(list_info)

        #peaks with only genes in forward strand
        elif gene_forward != "no_gene_plus" and gene_reverse == "no_gene_minus":
            #peaks where the forward gene doesn't have predicted function
            if gene_forward in dict_funct:
                list_info.insert(10,dict_funct[gene_forward])
                list_info.insert(13,'-')
                final_dict[peak] = tuple(list_info)
            else:
                list_info.insert(10,'No function detected')
                list_info.insert(13,'-')
                final_dict[peak] = tuple(list_info)



    return  final_dict


def removeChar(filename):
    file=os.path.splitext(os.path.basename(filename))[0] + "-annotationInfo.tsv"
    subprocess.call(['sed', '-i', 's/\"//g', file])

def main():

    parser = argparse.ArgumentParser(description='Script to analyse fold enrichment in ChipSeq peak files from macs2.')
    parser.add_argument(dest='peak_files', metavar='peakFiles', nargs='+', help='List of files where each line represents a peak file')
    parser.add_argument('--threshold', metavar= 'FLOAT', type=float, required = True, help='Fold enrichment threshold to filter peak files.')
    parser.add_argument('-s', '--sort', action='store_true', help = "Output new filtered files sorted by fold enrichment values above the threshold '-f'. Default: No output files will be written.")
    parser.add_argument('-gff', metavar='--gffFile',nargs=1,help = "Add genome GFF3 file to further analyse the regions of the genome with peaks above the threshold.")
    parser.add_argument('-funct', metavar='--functionalAnnotation', nargs=1, help= "Add blast like tab file with functional annotation of the genes predicted in the -gff file. [Support for m8 file from rapsearch2]")
    parser.add_argument('-col', metavar='--collumnToParse', nargs='+', type = int, help='Column/s number/s of the blastTAB file with functional annotation to parse. E.g. Rapsearch: Only column 2 [subject ID and description are in the same column]. \
     If more than one columns is desired, please set the numbers separated by space [e.g. "-col 2 3"].')
    args = parser.parse_args()

    if args.funct and not args.gff or args.funct and not args.col:
        logging.error("Error. '-funct' argument only available when '-gff' and '-col' is set.")
        exit(1)

    if args.col and not args.funct:
        logging.error("Error. '-col' argument only available when '-funct' is set.")
        exit(1)

    processFiles(args.peak_files,args.threshold, args.sort, args.gff,args.funct[0], args.col)

if __name__ == "__main__":
    main()