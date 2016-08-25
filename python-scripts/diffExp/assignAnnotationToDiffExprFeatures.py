__author__ = 'pedro'


import argparse
import logging
import sys
import os
import collections
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')

def processFeaturesIDs(listDiffExpressed,onlyIds ):
    dict = collections.OrderedDict() #dictionary where the feature ID will be the key and the rest a list of attributes dependding on the differential expression software used.
    commented_lines = []
    total_feature = 0
    logging.info("Processing " + listDiffExpressed + " file..")
    logging.info("Creating dictionary of differential expressed features..")
    with open(listDiffExpressed) as file:
        if onlyIds:
            for line in file:
                line.rstrip()
                if line.startswith("#"):
                    commented_lines.append(line)
                elif line.rstrip() in dict:
                    logging.error("Error Duplicate feature ID in " + listDiffExpressed + " file:\t" + line.rstrip())
                    exit(1)
                else:
                    dict[line.rstrip()] = []
                    total_feature += 1

        else:
            for line in file:
                line.rstrip()
                if line.startswith("#"):
                    commented_lines.append(line)
                else:
                    total_feature += 1
                    feature_fields = line.split()[1:]
                    measures = []
                    featureID = line.split()[0]
                    for field in feature_fields:
                        measures.append(field)
                    if featureID in dict:
                        logging.error("Warning - Duplicate feature ID in" + listDiffExpressed + " file:\t" + featureID )
                        dict[featureID].append(measures)

                    else:
                        dict[featureID] = [measures]
    file.close()
    return dict, total_feature, commented_lines

def processAnnotationFile(annotationFile, database_searched,dict_features, total_feature ,noDescription):
    annotated_features = 0
    logging.info("Processing annotation file ..")
    with open(annotationFile) as file:
        previous_query = ""
        if database_searched == "swissprot":
            for line in file:
                line.rstrip()
                if not line.startswith('#') :
                    query = line.split()[0]
                    if query in dict_features and query != previous_query:

                        annotated_features += 1
                        hit = line.split('\t')[1]
                        swissprot_id = hit.split("|")[1]
                        description = hit.split("|")[2]

                        #update dict
                        new_list = dict_features[query]
                        if len(new_list) != 0:
                             i=0
                             for measures in new_list:
                                 measures.extend([swissprot_id, description])
                                 new_list[i] = measures
                                 i+=1
                             dict_features[query] = new_list
                        else:
                            dict_features[query].append([swissprot_id, description])
                        previous_query = query


        elif database_searched == "ncbi-nr":
            for line in file:
                line.rstrip()
                if not line.startswith('#') :
                    query = line.split()[0]
                    if query in dict_features and query != previous_query:

                        annotated_features += 1
                        ncbi_id = line.split('\t')[1]
                        if noDescription:
                            description = "No description available"
                        else:
                            description = line.split('\t')[2]

                        #update dict
                        new_list = dict_features[query]
                        if len(new_list) != 0:
                             i=0
                             for measures in new_list:
                                 measures.extend([ncbi_id, description])
                                 new_list[i] = measures
                                 i+=1
                             dict_features[query] = new_list
                        else:
                            dict_features[query].append([ncbi_id, description])
                        previous_query = query


        elif database_searched == "eggnog":
            for line in file:
                line.rstrip()
                if not line.startswith('#') :
                    query = line.split()[2]
                    if query in dict_features and query != previous_query:

                        annotated_features += 1
                        eggnong_id = line.split("\t")[0].split(".")[1] #eggnog ortholog group in an hmmscan hits file
                        if not noDescription:

                            logging.error("When searching against eggnog (using hmmer), no description is available in the table output format (generated used "
                                          "tblout or domtblout or pfamtblout arguments on hmmscan. Please set the '-n' argument and rerun the script without descriptions.")
                            exit(1)

                        #update dict
                        new_list = dict_features[query]
                        if len(new_list) != 0:

                             i=0
                             for measures in new_list:
                                 measures.extend([eggnong_id])
                                 new_list[i] = measures
                                 i+=1
                             dict_features[query] = new_list
                        else:
                            dict_features[query].append([eggnong_id])
                        previous_query = query

    file.close()
    if annotated_features == 0:
        logging.error("No features found in the annotation file are present in the list of differential expressed genes. Please check if the feature IDs are concordant.")
        exit(1)


    percentage = round(float(annotated_features) / float(total_feature) * 100,2 )
    logging.info("Total number of differential expressed features in file\t" + str(total_feature))
    logging.info("Number of annotated features\t" + str(annotated_features))
    logging.info("Percentage annotated\t" + str(percentage) + '\n')
    return dict_features



def writeOutput(dict,commented_lines,outputFile,database):

    if os.path.exists(outputFile):
        os.remove(outputFile)
    with open(outputFile, 'w') as file:
        if len(commented_lines) > 0 and database == "eggnog":
            new_attributes_line = [commented_lines[-1].rstrip() + '\t' + database + '_id' + '\n']
            commented_lines = commented_lines[:-1] + new_attributes_line
            for line in commented_lines:
                file.write(line)
        elif len(commented_lines) > 0:
            new_attributes_line = [commented_lines[-1].rstrip() + '\t' + database + '_id' + '\t' + 'description' + '\n']
            commented_lines = commented_lines[:-1] + new_attributes_line
            for line in commented_lines:
                file.write(line)

        list_length_col = []
        max_col = 0
        for v in dict.values():
            for value in v:
                list_length_col.append(len(value))
            max_col = max(list_length_col)

        for k,v in dict.items():
            if v:
                for value in v:

                    if len(value) != max_col:
                        file.write(k + "\t" + "\t".join(value) + '\tNo annotation available\n')
                    else:
                        file.write(k + "\t" + "\t".join(value) + '\n')
            else:
                file.write(k + "\tNo annotation available\n")
    file.close()


def main():
    parser = argparse.ArgumentParser(description='Script to assign functional annotation to a list of differential expressed features based on a tab separated file '
                                                 'representing a similarity search. The features IDs must be concordant between files. The best hit approach will be selected.')
    parser.add_argument(dest='listDiffExpressed', metavar='listDiffExpressed', nargs=1, help='File where each line represents the significant features to process. Feature ID must be in the 1st columns.'
                                                ' All non-feature lines must start with a "#".')
    parser.add_argument(dest='annotationFile', metavar='annotationFile', nargs=1, help='Blast like tab separated file with functional annotations.')
    parser.add_argument(dest='output_file', metavar='output_file', nargs=1, help='File to write the output.')
    parser.add_argument(dest='database_searched', metavar='database_searched', nargs=1, type=str, choices=['swissprot','ncbi-nr', 'eggnog'],help='Database used for similarity searches. If ncbi-nr is selected, ncbi id must come on the 2nd column and its decription on the 3rd (See "-n" argument for exceptions). If eggnog choosen, it is supposed that '
        'hmmscan was run over hidden markov models representing the eggnog records. Available choices [swissprot, ncbi-nr, eggnog]')
    parser.add_argument(dest='software_used', metavar='software_used', nargs=1, type=str, choices=['rapsearch2', 'blastp','hmmer'],help='Tool used to run similarity searches. Available choices [rapsearch2, blastp]. When blastp used, it is assumed that the description comes in the 3rd columns of the annotation file. If not, set "-n" argument\
    to process only the IDs of the hits.')
    parser.add_argument('-id', '--idOnly', action='store_true', help='If set, listDiffExpressed represents file with only a feature ID per line, rather than a multi column file from EdgeR.')
    parser.add_argument('-n', '--noDescriptionBlastp', action='store_true', help='By default, blastp does not output database records description. If blastp was not run with special configuration to report this information in the\
    3rd column of the annotation file, please set this argumet')
    args = parser.parse_args()
	
   
    if "eggnog" in args.database_searched  and not "hmmer" in args.software_used:
        logging.error("Searches against eggnog database are only available using hmmer software. Please set 'hmmer' in the software used argument.")
        exit(1)
    elif "swissprot" in args.database_searched  or "ncbi-nr" in args.database_searched  and "hmmer" in args.software_used:
        logging.error("Searches against swissprot or any ncbi database (except for CDD domain database) can not be performed with hmmer software. Please change"
                      " the software used argument to a different value.")
        exit(1)


    dict_features, total_features, commented_lines= processFeaturesIDs(args.listDiffExpressed[0], args.idOnly)
    dict_updated = processAnnotationFile(args.annotationFile[0], args.database_searched[0], dict_features,total_features, args.noDescriptionBlastp)
    writeOutput(dict_updated,commented_lines, args.output_file[0], args.database_searched[0])

if __name__ == "__main__":
    main()
