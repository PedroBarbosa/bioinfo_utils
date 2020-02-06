__author__ = 'pedro'

import argparse
import os
import logging
import sys
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')

def processMAC2xls(xlsPeakFile, basename):
    logging.info("Reading and creating dictionary for %s peaks file.." % xlsPeakFile)
    dict = {}

    output_xls = basename + '_MAC2_IDR_peaks.xls'
    if os.path.exists(output_xls):
        os.remove(output_xls)


    with open(output_xls, 'w') as file_out:

        with open(xlsPeakFile, 'r') as file_in:
            for line in file_in:
                if line.startswith("#") or "fold_enrichment" in line:
                    file_out.write(line)
                else:
                    line_features = line.split()
                    key = ':'.join(line_features[0:3])
                    dict[key] = tuple(line_features[3:])

    file_in.close()
    file_out.close()
    return dict, output_xls

def processIDRnarrowPeak(idrFile, dict_macs2_xls, outputXLS_file, basename):

    listIDR_peaks = []
    output_ODR  = basename + '_IDRformat.narrowPeak'
    if os.path.exists(output_ODR):
        os.remove(output_ODR)

    logging.info("Writing new XLS file with only the peaks that passed the IDR threshold and new IDR narrowPeak file with peak names and fold enrichment values added.")
    with open(outputXLS_file, 'a') as file_xls_out:
        with open(output_ODR, 'w') as file_idr_out:
            with open(idrFile, 'r') as file_idr_in:

                for line in file_idr_in:
                    line_features = line.strip().split()
                    key = [line_features[0], int(line_features[1]) + 1, line_features[2]]
                    new_list = [str(k) for k in key]
                    key_merged = ':'.join(new_list)


                    if key_merged in dict_macs2_xls:
                        #write to new xls only peaks that passed IDR thresholds
                        peak_coord = key_merged.split(":")
                        file_xls_out.write(peak_coord[0] + '\t' + peak_coord[1] + '\t' + peak_coord[2] + '\t' + '\t'.join(dict_macs2_xls[key_merged]) + '\n')

                        #new IDR peak narrowPeak file with peak names and foldEnrichment values added to the last columns
                        peak_name = dict_macs2_xls[key_merged][-1]
                        peak_foldEnrichment = dict_macs2_xls[key_merged][4]
                        listIDR_peaks.append(peak_name)

                        if line_features[3] == '.':
                            line_features[3] = peak_name
                        line_features.append(peak_foldEnrichment)
                        file_idr_out.write('\t'.join(line_features) + '\n')

                    else:
                        logging.warning("Scaffold or peak coordinates of the %s IDR peak were not found in MACS2 xls file." % key_merged)


    file_idr_out.close()
    file_xls_out.close()
    file_idr_in.close()
    return listIDR_peaks


def processOriginalNarrowPeakFromMACS2(originalNarrowPeak, listIDRpeaks, basename):
    logging.info("Processing original MACS2 narrowPeak file.")
    dict = {}
    output_originalNarrow = basename + '_MAC2_IDR_peaks.narrowPeak'
    if os.path.exists(output_originalNarrow):
        os.remove(output_originalNarrow)


    with open(output_originalNarrow, 'w') as file_narrowPeak_out:
        with open(originalNarrowPeak, 'r') as file_narrowPeak_in:
            for line in file_narrowPeak_in:
                if line.startswith("track type"):
                    file_narrowPeak_out.write(line.rstrip() + '\n')
                else:
                    line_attributes = line.strip().split()
                    peak_name =line_attributes.pop(3)
                    dict[peak_name] = line_attributes

        file_narrowPeak_in.close()

        logging.info("Writing new narrowPeak file with the peaks that passed the IDR threshold.")
        for peak in listIDRpeaks:
            if peak in dict:
                attr = dict[peak]
                attr.insert(3,peak)
                file_narrowPeak_out.write('\t'.join(attr) + '\n')
            else:
                logging.warning("IDR Peak %s not present in original MACS2 narrowPeak file." % peak)

    file_narrowPeak_out.close()




def main():

    parser = argparse.ArgumentParser(description='Script to process peaks that passed the IDR threshold (version2) based on the xls outputted by macs2. Starts from the\
     assumption that fold enrichment values and peak names in the IDR edited narrowPeak format are missing. It outputs the IDR narroPeak format with peak names in 4th columns and adds\
                                                 fold enrichment values to the last column. Furthermore, a filtered xls and optionally, a narrowPeak file will be written\
                                                 displaying only the peaks that passed the IDR threshold.')
    parser.add_argument(dest='IDR_peakFile', metavar='idrPeakFile', nargs=1, help='Extended narrowPeak file outputted by IDR v2.0.')
    parser.add_argument(dest='macs2_xls_file', metavar='xlsPeaksFile', nargs=1, help='XLS file produced by macs2, should be analogous to the narrowPeak file that served\
    as the oracle file when IDR was ran.')
    parser.add_argument('--originalNarrowPeak', metavar= 'narrowPeak', help='Original narrowPeak file (from the same run of "xlsPeaksFile"). Use this if you wants to ouput \
    the original narrowPeak file with the peaks that passed the IDR threshold.')
    parser.add_argument('-o', metavar='--outputBasename', required=True, help = "[Required] Basename to write the output files.")
    args = parser.parse_args()

    dict_xls, outputFile_xls = processMAC2xls(args.macs2_xls_file[0], args.o)
    listIDRpeaks = processIDRnarrowPeak(args.IDR_peakFile[0], dict_xls,outputFile_xls, args.o)
    if args.originalNarrowPeak:
        processOriginalNarrowPeakFromMACS2(args.originalNarrowPeak, listIDRpeaks, args.o)



if __name__ == "__main__":
    main()