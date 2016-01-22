__author__ = 'pedro'

import argparse
import sys
import os
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
import collections


def processTxtTab(inputFiles):
    logging.info("Processing tab separated txt files..")
    mydict = {}
    diffexpressed_features_per_comparison = collections.OrderedDict()
    for filename in inputFiles:

        with open(filename) as file:
            comparison = os.path.splitext(os.path.basename(filename))[0]
            numb_features = 0
            for line in file:
                if not line.startswith('#'):
                    line.rstrip()
                    numb_features += 1
                    featureID = line.split('\t')[0]

                    if featureID not in mydict:
                        overlapping_comparisons_perFeature = []
                        overlapping_comparisons_perFeature.append(comparison)
                        mydict[featureID] = overlapping_comparisons_perFeature

                    else:
                        #update comparison file
                        common_comparisons = mydict[featureID]
                        common_comparisons.append(comparison)
                        mydict[featureID] = common_comparisons

            diffexpressed_features_per_comparison[comparison] = numb_features

    return mydict,diffexpressed_features_per_comparison



def processFromIDs(inputFiles):
    logging.info("Processing feature identifiers files.. ")
    mydict = {}
    diffexpressed_features_per_comparison = collections.OrderedDict()
    for filename in inputFiles:

        with open(filename) as file:
            comparison = os.path.splitext(os.path.basename(filename))[0]
            numb_features = 0
            for featureID in file:
                featureID.rstrip()
                numb_features += 1
                if featureID not in mydict:
                    overlapping_comparisons_perFeature = []
                    overlapping_comparisons_perFeature.append(comparison)
                    mydict[featureID] = overlapping_comparisons_perFeature

                else:
                    #update comparison file
                    common_comparisons = mydict[featureID]
                    common_comparisons.append(comparison)
                    mydict[featureID] = common_comparisons

            diffexpressed_features_per_comparison[comparison] = numb_features

    return mydict,diffexpressed_features_per_comparison







def writeOutput(mydict,features_per_comparison):#,outputFile):
    final_dict =  collections.OrderedDict()
    ocurrences_dict = collections.OrderedDict()

#    if os.path.exists(outputFile):
#        os.remove(outputFile)

    logging.info("Processing intersections and printing output ..")
#    with open(outputFile, "w") as file:

    print('#Comparison' + '\t' + '#Number of differential features' + '\n')
    for comparison,number_features in features_per_comparison.iteritems():
        print(comparison + '\t' + str(number_features) + '\n')
    print('\n\n')
    print('#Total number of unique features with differential expression in any comparison' + '\t' + str(len(mydict)))
    print('\n\n\n')


    #Create dict relating number of times the comparison have differential expression per feature
    for i in range(1,len(features_per_comparison.keys()) + 1,1):
        ocurrences_dict[i] = 0


    string_binary =""
    final_dict['#FeatureID'] = ['Comparisons with diffExp', '\t'.join(features_per_comparison.keys())]
    for featureID, comparisons in mydict.iteritems():
        for comparison in features_per_comparison.keys():
            if comparison in comparisons:
                string_binary = string_binary + 'yes' + '\t'
            else:
                string_binary = string_binary + 'no' + '\t'


        ocurrences_dict[len(mydict[featureID])] += 1
        final_dict[featureID] = [str(len(mydict[featureID])), string_binary ]
        string_binary = ""


    print('#Table displaying number of times each feature has N number of comparisons with differential expression.\n')
    for k,v in ocurrences_dict.iteritems():
        print(str(k) + ' comparisons' + '\t' + str(v) + '\n')


    file.write('\n\n\n#Table displaying the comparisons in which the features have differential expression.\n')
    for k,v in final_dict.iteritems():
        print(k + '\t' + v[0] + '\t' + v[1].rstrip() + '\n')

#    file.close()
    return final_dict, ocurrences_dict



def main():

    parser = argparse.ArgumentParser(description='Script to check the interception of the differential expressed features present in different pairwise tests.')
    parser.add_argument(dest='input_files', metavar='diffExp_files', nargs='+', help='List of files where each line represents the significant features to process. Feature ID must be in the 1st column.'
                                                    ' All non-feature lines must start with a "#" (minimum 2 files).')
#    parser.add_argument('-o', "--output", required =True, help='File to write the output.')
    parser.add_argument('-l', '--list', action='store_true', help='Process feature identifiers (one per line) rather than txt tab separated output files.')
    args = parser.parse_args()


###############################################################################
########################## COMMAND LINE PARSING ###############################
##############################################################################

    if len(args.input_files) < 2 :
        logging.fatal('Error: %s\n' % 'You should specify at least two files to compare.')

    elif args.list:
        mydict,diffexpressed_features_per_comparison=processFromIDs(args.input_files)
        writeOutput(mydict,diffexpressed_features_per_comparison)#, args.output)
        #intersection(dict_unique,dict_repeated)

    else:
        mydict,diffexpressed_features_per_comparison = processTxtTab(args.input_files)
        writeOutput(mydict,diffexpressed_features_per_comparison)#, args.output)
        #intersection(dict_unique,dict_repeated)


if __name__ == "__main__":
    main()