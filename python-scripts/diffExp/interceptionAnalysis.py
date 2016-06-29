__author__ = 'pedro'

import argparse
import sys
import os
import logging
logging.basicConfig(stream=sys.stdout, level=logging.ERROR, format='%(asctime)s %(message)s')
import collections


def processWithlogFC(inputFiles):
    logging.info("Processing tab separated txt files..")
    mydict = {}
    diffexpressed_features_per_comparison = collections.OrderedDict()
    for filename in inputFiles:

        with open(filename) as file:
            comparison = os.path.splitext(os.path.basename(filename))[0]
            numb_features = 0
            for line in file:
                line.replace(" ","\t")
                if not line.startswith('#'):
                    line = line.rstrip()
                    numb_features += 1
                    featureID = line.split('\t')[0]
                    logFC = line.split('\t')[1]
                    if featureID not in mydict:
                        overlapping_comparisons_perFeature = []
                        if float(logFC) > 0:
                            tuple = (comparison,logFC,'up')
                            overlapping_comparisons_perFeature.append(tuple)
                            mydict[featureID] = overlapping_comparisons_perFeature

                        else:
                            tuple = (comparison,logFC,'down')
                            overlapping_comparisons_perFeature.append(tuple)
                            mydict[featureID] = overlapping_comparisons_perFeature

                    else:
                        #update comparison file
                        common_comparisons = mydict[featureID]
                        if float(logFC) > 0:
                            tuple = (comparison,logFC,'up')
                            common_comparisons.append(tuple)
                            mydict[featureID] = common_comparisons
                        else:
                            tuple = (comparison,logFC,'down')
                            common_comparisons.append(tuple)
                            mydict[featureID] = common_comparisons

            diffexpressed_features_per_comparison[comparison] = numb_features

    return mydict,diffexpressed_features_per_comparison


def processTxtTab(inputFiles):
    logging.info("Processing tab separated txt files..")
    mydict = {}
    diffexpressed_features_per_comparison = collections.OrderedDict()
    for filename in inputFiles:

        with open(filename) as file:
            comparison = os.path.splitext(os.path.basename(filename))[0]
            numb_features = 0
            for line in file:
                line.replace(" ","\t")
                if not line.startswith('#'):
                    line = line.rstrip()
                    numb_features += 1
                    featureID = line.split('\t')[0]
                    print(featureID)
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
                featureID = featureID.rstrip()
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




def writeOutputWithLogFC(mydict,features_per_comparison):#,outputFile):
    final_dict =  collections.OrderedDict()
    ocurrences_dict = collections.OrderedDict()

#    if os.path.exists(outputFile):
#        os.remove(outputFile)

    logging.info("Processing intersections and printing output ..")
#    with open(outputFile, "w") as file:

    print('#Comparison' + '\t' + '#Number of differential features')
    for comparison,number_features in iter(features_per_comparison.tems()):
        print(comparison + '\t' + str(number_features))
    print('\n')
    print('#Total number of unique features with differential expression in any comparison' + '\t' + str(len(mydict)))
    print('\n\n')


    #Create dict relating number of times the comparison have differential expression per feature
    for i in range(1,len(features_per_comparison.keys()) + 1,1):
        ocurrences_dict[i] = 0


    string_row =""
    string_columns = ""
    for comparison in features_per_comparison.keys():
        string_columns += comparison + '\t' + 'logFC [' + comparison + ']\t'


    final_dict['#FeatureID'] = ['Comparisons with diffExp', string_columns ]
    for featureID, comparisons in iter(mydict.tems()):

        for comparison in features_per_comparison.keys():
            l = [item for item in comparisons if comparison in item]
            if not l:
                string_row = string_row + 'noDiff' + '\t' + '-' + '\t'
            else:
                tuple = [i[1] for i in l]
                string_row = string_row + ''.join([i[2] for i in l]) + '\t' + ''.join([i[1] for i in l]) + '\t'



        ocurrences_dict[len(mydict[featureID])] += 1
        final_dict[featureID] = [str(len(mydict[featureID])), string_row ]
        string_row = ""


    print('#Table displaying number of times each feature has N number of comparisons with differential expression.')
    for k,v in iter(ocurrences_dict.items()):
        print(str(k) + ' comparisons' + '\t' + str(v))

#    file.close()

    return final_dict, ocurrences_dict


def writeOutput(mydict,features_per_comparison):#,outputFile):
    final_dict =  collections.OrderedDict()
    ocurrences_dict = collections.OrderedDict()

#    if os.path.exists(outputFile):
#        os.remove(outputFile)

    logging.info("Processing intersections and printing output ..")
#    with open(outputFile, "w") as file:

    print('#Comparison' + '\t' + '#Number of differential features')
    for comparison,number_features in iter(features_per_comparison.items()):
        print(comparison + '\t' + str(number_features))
    print('\n')
    print('#Total number of unique features with differential expression in any comparison' + '\t' + str(len(mydict)))
    print('\n\n')


    #Create dict relating number of times the comparison have differential expression per feature
    for i in range(1,len(features_per_comparison.keys()) + 1,1):
        ocurrences_dict[i] = 0


    string_binary =""
    final_dict['#FeatureID'] = ['Comparisons with diffExp', '\t'.join(features_per_comparison.keys())]
    for featureID, comparisons in iter(mydict.items()):
        for comparison in features_per_comparison.keys():
            if comparison in comparisons:
                string_binary = string_binary + 'yes' + '\t'
            else:
                string_binary = string_binary + 'no' + '\t'


        ocurrences_dict[len(mydict[featureID])] += 1
        final_dict[featureID] = [str(len(mydict[featureID])), string_binary ]
        string_binary = ""


    print('#Table displaying number of times each feature has N number of comparisons with differential expression.')
    for k,v in iter(ocurrences_dict.items()):
        print(str(k) + ' comparisons' + '\t' + str(v))


#    file.close()
    return final_dict, ocurrences_dict

def printFinalWithAnnotation(final_dict,annotationFile):

    output_dict = processFromBlastTab(final_dict, annotationFile)

    print('\n\n#Table displaying the comparisons in which the features have differential expression.')
    header = output_dict['#FeatureID']
    header.extend(['Annotation ID','Description'])
    output_dict['#FeatureID'] = header

    for k,v in iter(output_dict.items()):
        if len(v) > 2:
            print(k + '\t' + v[0] + '\t' + v[1].rstrip() + '\t' + v[2] + '\t' + v[3])
        else:
            print(k + '\t' + v[0] + '\t' + v[1].rstrip())





def printFinal(final_dict):
    print('\n\n#Table displaying the comparisons in which the features have differential expression.')
    for k,v in iter(final_dict.items()):
        print(k + '\t' + v[0] + '\t' + v[1].rstrip())





def processFromBlastTab(dict_final, annotationFile):
    annotated_features = 0
    previous_query = ""
    swissprot = False
    ncbi_nr = False
    database_id,description = "",""
    with open(annotationFile[0], 'r') as annotFile:
        for line in annotFile:
            line.rstrip()
            if not line.startswith('#') :
                query = line.split()[0]
                if query in dict_final:

                    if query == previous_query:
                        previous_query = query
                    else:
                        annotated_features += 1
                        hit = line.split('\t')[1]
                        if swissprot:
                            database_id = hit.split("|")[1]
                            description = hit.split("|")[2]
                        elif ncbi_nr:
                            protein_id = hit.split("|")[1]
                            refseq_id = hit.split("|")[3]
                            database_id = "NCBI_id:" + protein_id + ",RefSeq_id:" + refseq_id
                            description = hit.split("|")[4]
                            description = description.split("]")[0]
                            description = description + "]"

                        #update dict
                        new_list = dict_final[query]
                        new_list.extend([database_id, description])
                        dict_final[query] = new_list
                        previous_query = query
            elif 'nr' in line:
                ncbi_nr = True

            elif 'swissprot' in line:
                swissprot = True

    if annotated_features == 0:
        logging.error("No features found in the annotation file are present in the list of differential expressed genes. Please check if the feature IDs are concordant.")
        exit(1)


    percentage = round(float(annotated_features) / float(len(dict_final) - 1) * 100,2 )
    print ("\n\nNumber of annotated features\t" + str(annotated_features))
    print ("Percentage annotated\t" + str(percentage) + '\n')
    return dict_final



def main():

    parser = argparse.ArgumentParser(description='Script to check the interception of the differential expressed features present in different pairwise tests.')
    parser.add_argument(dest='input_files', metavar='diffExp_files', nargs='+', help='List of files where each line represents the significant features to process. Feature ID must be in the 1st column.'
                                                    ' All non-feature lines must start with a "#" (minimum 2 files).')
#    parser.add_argument('-o', "--output", required =True, help='File to write the output.')
    parser.add_argument('-l', '--list', action='store_true', help='Process feature identifiers (one per line) rather than txt tab separated output files.')
    parser.add_argument('-a', '--all', action='store_true', help='Include log fold changes in the ouptut. Tab separated files are required. By default (edgeR), logFC values appear in the second column.')
    parser.add_argument('-f', dest='annotationFile', help='Add annotations to the output. Please add a file with the annotations in the blastTAB format.')
    args = parser.parse_args()


###############################################################################
########################## COMMAND LINE PARSING ###############################
###############################################################################

    if len(args.input_files) < 2 :
        logging.fatal('Error: %s\n' % 'You should specify at least two files to compare.')

    elif args.list and args.all:
        logging.fatal('Error: %s\n' % 'If --all provided, txt tab separated files are required. Please remove --list parameter.')

    elif args.list:
        mydict,diffexpressed_features_per_comparison=processFromIDs(args.input_files)
        writeOutput(mydict,diffexpressed_features_per_comparison)#, args.output)


    elif args.all:
        mydict,diffexpressed_features_per_comparison = processWithlogFC(args.input_files)
        final_dict, ocurrencesDict = writeOutputWithLogFC(mydict,diffexpressed_features_per_comparison)#, args.output)
        if args.annotationFile:
            printFinalWithAnnotation(final_dict,args.annotationFile)
        else:
            printFinal(final_dict)

    else:
        mydict,diffexpressed_features_per_comparison = processTxtTab(args.input_files)
        final_dict, ocurrencesDict = writeOutput(mydict,diffexpressed_features_per_comparison)#, args.output)
        if args.annotationFile:
            printFinalWithAnnotation(final_dict,args.annotationFile)
        else:
            printFinal(final_dict)


if __name__ == "__main__":
    main()