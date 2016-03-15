__author__ = 'pedro'

import argparse
import logging
import sys
import os
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
from collections import defaultdict


def getOverlap(list_dict,list_files, list_number_peaks):
    #create list of sets
    set_list = []
    for dict in list_dict:
        s1 = set(dict.keys())
        set_list.append(s1)

    #set of scaffold_id in common among all comparisons
    interception = set.intersection(*set_list)

    for i in range(0,len(list_files), 1):
        print "Number of peaks in file %s:\t%i" % (os.path.basename(list_files[i]),list_number_peaks[i][0])
        print "Number of different scaffolds with peaks in file %s:\t%i\n" % (os.path.basename(list_files[i]),list_number_peaks[i][1])

    print "\nMaximum number of peak overlaps:\t%s" % min(list_number_peaks, key = lambda m: m[0])[0]
    print "Number of possible peaks with overlaps (peaks in scaffolds shared across all samples):\t%i\n" % len(interception)
    overlaps = 0
    #coordinates analysis
    for scaffold_id in interception:
        ref_dict = list_dict[0]
        coord = ref_dict[scaffold_id]
        for i in range(1,len(list_dict),1):
            temp_dict = list_dict[i]
            coord2 = temp_dict[scaffold_id]

            for peak in coord:
                for peak2 in coord2:

                    if peak[2] >= peak2[1] and peak[1] <= peak2[1]: #if end of peak1 is bigger than beginning of peak2
                        overlaps += 1
                        file_ref = os.path.basename(list_files[0])
                        file_com = os.path.basename(list_files[i])
                        logging.info("'%s' and '%s' peaks in '%s' and '%s' files do overlap somehow." % (peak[0], peak2[0], file_ref, file_com))
                    elif peak[1] <= peak2[2] and peak[1] >= peak2[1]: #if beggining of peak1 is smaller thanafter of peak2
                        overlaps += 1
                        file_ref = os.path.basename(list_files[0])
                        file_com = os.path.basename(list_files[i])
                        logging.info("'%s' and '%s' peaks in '%s' and '%s' files do overlap somehow." % (peak[0], peak2[0], file_ref, file_com))
                    elif peak[1] >= peak2[1] and peak[1] <= peak2[2]: #if peak 1 total within
                        overlaps += 1
                        file_ref = os.path.basename(list_files[0])
                        file_com = os.path.basename(list_files[i])
                        logging.info("'%s' and '%s' peaks in '%s' and '%s' files do overlap somehow." % (peak[0], peak2[0], file_ref, file_com))
                    elif peak2[1] >= peak[1] and peak2[2] <= peak[2]: #if peak 2 total within
                        overlaps += 1
                        file_ref = os.path.basename(list_files[0])
                        file_com = os.path.basename(list_files[i])
                        logging.info("'%s' and '%s' peaks in '%s' and '%s' files do overlap somehow." % (peak[0], peak2[0], file_ref, file_com))

    print "Number of overlapped peaks found:\t%i" % overlaps



def processFiles(peak_files):


    list_files = []
    list_dicts = []
    list_number_peaks = []
    for filename in peak_files:
        peaks = 0
        peaks_different_scaff = 0
        dict = defaultdict(list)
        list_files.append(filename)
        with open(filename,'r') as file:

            for line in file:
                if not line.startswith("#"):
                    peaks += 1
                    line = line.rstrip()
                    fields = line.split("\t")
                    if len(fields) != 4:
                        logging.error("Input file %s does not comply with the requirements: 4 columns. Call help message on the script." % filename)
                    elif fields[1] not in dict:
                        peaks_different_scaff += 1
                        dict[fields[1]] = [(fields[0],fields[2], fields[3])]
                    else:
                        dict[fields[1]].append((fields[0],fields[2], fields[3]))

            list_dicts.append(dict)
            list_number_peaks.append((peaks,peaks_different_scaff))

        file.close()

    return list_dicts,list_files, list_number_peaks


def main():

    parser = argparse.ArgumentParser(description='Script to analyse common peaks  enrichment in ChipSeq peak files from macs2.')
    parser.add_argument(dest='peak_files', metavar='peakFiles', nargs='+', help='List of files where each line represents a peak. Must come in a tab separated 4 columns format: peak_id, scaffold_id, start_reference, stop_reference.')
    args = parser.parse_args()

    if len(args.peak_files) < 2:
        logging.error("Error.Please add at least two files to compare!")
    else:
        list_dics, list_files, list_number_peaks = processFiles(args.peak_files)
        getOverlap(list_dics, list_files,list_number_peaks)

if __name__ == "__main__":
    main()