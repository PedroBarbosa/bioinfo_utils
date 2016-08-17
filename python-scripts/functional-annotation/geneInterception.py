__author__ = 'pedro'

import argparse
import sys
import os
import csv
from itertools import combinations
from collections import defaultdict
from collections import OrderedDict


def processFromBlastTab(inputFiles,bestHit):
    mydict_unique = defaultdict(list) #dictionary of all unique IDs obtained annotation file
    mydict_repeated = defaultdict(list) #dictionary of the repeated IDs per annotation file
    for filename in inputFiles:
        list_name = filename
        mydict_unique[list_name] = []
        mydict_repeated[list_name] = []
        annotated_reads = 0


        with open(filename) as file:

            previous_query = ""
            print("Processing file " + filename + "..")
            if bestHit:

                for line in file:
                    if not line.startswith('#') and len(line.split()) == 1:
                        print("%s file represents gene identifiers, please set '-g' argument accordingly. Exiting.." % filename)
                        exit(2)
                    if not line.startswith('#') and line.rstrip():

                        query = line.split()[0]
                        if query != previous_query:
                            annotated_reads += 1
                            hit = line.split()[1]
                            values = mydict_unique[list_name]
                            mydict_unique[list_name].append(hit) if hit not in values else mydict_repeated[list_name].append(hit)

                        previous_query = query


            else:
                for line in file:
                    if not line.startswith("#") and len(line.split()) == 1:
                        print("%s file represents gene identifiers, please set '-g' argument accordingly. Exiting.." % filename)
                        exit(2)

                    if not line.startswith('#'):

                        query = line.split()[0]
                        if query != previous_query:
                            annotated_reads += 1

                        hit = line.split()[1]
                        values = mydict_unique[list_name]
                        mydict_unique[list_name].append(hit) if hit not in values else mydict_repeated[list_name].append(hit)
                        previous_query = query

        print("%s annotated genes/reads in file %s!" % (annotated_reads,filename))
    return annotated_reads,mydict_unique,mydict_repeated


def processFromHmmer(inputFiles,bestHit,eggNOG):
    mydict_unique = defaultdict(list) #dictionary of all unique IDs obtained annotation file
    mydict_repeated = defaultdict(list) #dictionary of the repeated IDs per annotation file
    for filename in inputFiles:
        list_name = filename
        mydict_unique[list_name] = []
        mydict_repeated[list_name] = []
        annotated_reads = 0
        with open(filename) as file:
            previous_query = ""
            print("Processing file " + filename + "..")

            if bestHit:
                for line in file:
                    if not line.startswith('#') and line.rstrip():
                        query = line.split()[2]
                        if query != previous_query:

                            annotated_reads += 1
                            if eggNOG:
                                hit = line.split("\t")[0].split(".")[1] #eggnog ortholog group in an hmmscan hits file
                            else:
                                hit = line.split()[0]

                            values = mydict_unique[list_name]
                            mydict_unique[list_name].append(hit) if hit not in values else mydict_repeated[list_name].append(hit)
                        previous_query = query


            else:
                for line in file:
                    if not line.startswith('#'):
                        query = line.split()[2]
                        if query != previous_query:
                            annotated_reads += 1

                        if eggNOG:
                            hit = line.split()[0].split(".")[1] #eggnog ortholog group in an hmmscan hits file
                        else:
                            hit = line.split()[0]

                        values = mydict_unique[list_name]
                        mydict_unique[list_name].append(hit) if hit not in values else mydict_repeated[list_name].append(hit)
                        previous_query = query


        print("%s annotated genes/reads in file %s!" % (annotated_reads,filename))
    return annotated_reads,mydict_unique,mydict_repeated


def processFromIDs(inputFiles):
    #createDicts
    mydict_unique = defaultdict(list)
    mydict_repeated = defaultdict(list)

    for filename in inputFiles:

        with open(filename) as file:
            list_name = filename
            mydict_unique[list_name] = []
            mydict_repeated[list_name] = []
            for gi in file:
                values = mydict_unique[list_name]
                GI = gi.rstrip()
                mydict_unique[list_name].append(GI) if GI not in values else mydict_repeated[list_name].append(GI)


    return mydict_unique,mydict_repeated

def intersection(dict_uniq,dict_repeat,outputBasename, descriptionFile):
    #Get interception list
    dic_annDescription = {}
    print("\nCalculating interception..")
    print("Generating stats..\n")

    interception = set()

    for v in iter(dict_uniq.values()):
        interception.update(v)

    if descriptionFile:
        with open(descriptionFile, 'r') as annotFile:

            for line in annotFile:
                if len(line.split("\t")) == 2:
                    id = line.split("\t")[0].rstrip()
                    if id in interception:
                        dic_annDescription[id] = line.split("\t")[1].rstrip()
                else:
                    print("Annotations description file is not tab separated or doesn't have 2 columns. Please check the format of the file or just remove '-a' argument. Exiting..")
                    exit(2)
            if len(dic_annDescription) == 0:
                print("No gene ID was concordant between the annotations files and the description (-a) file. Descriptions will not be added.")

    max_possible = 100000000000
    smallest_set = ""
    list_of_sets, list_of_k = [],[]
    i = 0
    total_ids = len(interception)
    for k,v in iter(dict_uniq.items()):

        interception = interception & set(v)
        if len(v) < max_possible:
            max_possible = len(v)
            smallest_set = k.split("/")[-1]

        list_of_k.append(k)
        list_of_sets.append((i,set(v)))
        i += 1

    print("Total number of existing IDs across all annotations:\t%s" % total_ids)
    print ("Maximum number of possible overlaps (Size of the smallest set of IDs in all files):\t%s (%s)" % (max_possible,smallest_set))
    print ("Number of annotations shared among all:\t%s (%s)" % (len(interception),round(float(len(interception)) / float(max_possible) * 100,2)) + '\n')

    #Percentage
    for k,v in iter(dict_uniq.items()):
        percentage = round(float(len(interception)) / float((len(v))) * 100,2)
        print("Percentage of %s shared across all annotations:\t%s" % (k.split("/")[-1], percentage))

    #Pairwise comparison
    print("\n")
    for i in combinations(list_of_sets, 2):
        inter = i[0][1] & i[1][1]
        print("Number of annotations shared between %s and %s:\t%s" % (list_of_k[i[0][0]].split("/")[-1], list_of_k[i[1][0]].split("/")[-1], str(len(inter))))


    ##Unique analysis
    newDict=defaultdict(set)
    print("\n")
    print("#################Single copy genes analysis ################\n#########################################################")
    for k,v in iter(dict_uniq.items()):
        [newDict[id].add(k.split("/")[-1]) for id in v]



    print("Writing unique IDs output file to %s directory.." % os.path.basename(outputBasename))
    with open(outputBasename + "_uniqueIDs.txt", "w") as csvfile:
        writer = csv.writer(csvfile,dialect=csv.excel_tab)
        writer.writerow(("#List of unique IDs to each annotation",''))

        dict_IDs_1annotation = defaultdict(list)
        for id,samples in iter(newDict.items()):
            if len(samples) == 1: #check if only one annotation has this id
                sample = next(iter(samples))
                dict_IDs_1annotation[sample].append(id)


        for k in sorted(dict_IDs_1annotation, key=lambda k: len(dict_IDs_1annotation[k]), reverse=True):
            [writer.writerow((k,id)) for id in dict_IDs_1annotation[k]]
    csvfile.close()



    print("Writing intersections to %s directory..\n\n" % os.path.basename(outputBasename))
    with open(outputBasename + "_interseptionTable.txt", "w") as csvfile:
        #writer = csv.writer(csvfile,dialect=csv.excel_tab,escapechar='',doublequote=False, quoting=csv.QUOTE_NONE)
        i = 0
        listKeys = []
        dict_withIndex = OrderedDict()
        for k in sorted(dict_uniq.keys()):
            key = k.split("/")[-1]
            listKeys.append(key)
            dict_withIndex[key] = i
            i += 1

        if len(dic_annDescription) > 0:
            csvfile.write("#List of unique IDs to each annotation" + "\t" + "\t".join(listKeys) + "\tFeature description\n")
            #writer.writerow(("#List of unique IDs to each annotation",'\t'.join(listKeys), "Feature description"))
        else:
            csvfile.write("#List of unique IDs to each annotation" + "\t" + "\t".join(listKeys) + "\n")
            #writer.writerow(("#List of unique IDs to each annotation",'\t'.join(listKeys)))


        for k in sorted(newDict, key=lambda k: len(newDict[k]), reverse=True):

            final_list = ['no'] * len(listKeys)
            for annotations in newDict[k]:
                index = dict_withIndex[annotations]
                final_list[index] = 'yes'

            if len(dict_IDs_1annotation) > 0 and k in dic_annDescription.keys():
                csvfile.write(k + "\t" + "\t".join(final_list) + "\t" + dic_annDescription[k] + "\n")
                #writer.writerow((k + '\t'.join(final_list),dic_annDescription[k]))


            else:
                csvfile.write(k + "\t" + "\t".join(final_list)+ "\n")
                #writer.writerow((k,'\t'.join(final_list)))

        #removeChar(outputBasename + "_interseptionTable.txt")

    csvfile.close()

    #Repeated genes analysis
    newDict_rep = defaultdict(list)
    print("\n")
    print("####################Repeated genes analysis ################\n################################################################")
    #populate dict
    all_repeated_ids = set()
    for v in iter(dict_repeat.values()):
        all_repeated_ids.update(v)

    for k in iter(dict_repeat.keys()):
        for id in all_repeated_ids:
            newDict_rep[id].append([k.split("/")[-1],0])

    print("Total number of genes with more than one copy in any of the genomes:\t%s\n" % len(all_repeated_ids))

    for k,v in iter(dict_repeat.items()):

        print("Number of genes with more than one copy in the genome for %s sample:\t%s" % (k.split("/")[-1],len(set(v))))

        for repeated in v:
            for element in newDict_rep[repeated]:
                ind = newDict_rep[repeated].index(element)
                if element[0] == k.split("/")[-1] and element[1] == 0:
                    newDict_rep[repeated][ind] = [element[0],2]
                elif element[0] == k.split("/")[-1]:
                    newDict_rep[repeated][ind] = [element[0], element[1] + 1]

            #print(a)
           # if newDict_rep[repeated][k.split("/")[-1]] == 0:
           #     newDict_rep[repeated][k.split("/")[-1]] += 2
           # else:
           #     newDict_rep[repeated][k.split("/")[-1]] += 1


    print("\nFrequency of the > 1x copy genes across samples:")
    if bool(newDict_rep):

        for key,value in sorted(newDict_rep.items(), key=lambda e: e[1][0:len(e)], reverse=True): # sorted(newDict_rep,lambda v: newDict_rep[v][0]):#, reverse=True):
            #print(value[0:len(value)])
            print(key,value)
        #for id in sorted(list(newDict_rep.values()),key=lambda (k,p) :(newDict_rep[p][k]),reverse=True):#,key=lambda k: k[id],reverse=True):# key=lambda x: x[0], reverse =True):#lambda x: newDict_rep[x][0], reverse=True):
        #    print(id, newDict_rep[id])
        #    break
    else:
        print("No genes found on this condition.")


parser = argparse.ArgumentParser(description='Script to check the interception of the genes present in different genome annotations (default blastTAB format).')
parser.add_argument(dest='input_files', metavar='annotated_files', nargs='+', help='Annotation files to be analyzed (minimum 2).')
parser.add_argument('-o', metavar='outputBasename', required = True, help='Output basename for which the output files will be written.')
parser.add_argument('-g', '--GIlist', action='store_true', help='Process gene identifiers (one per line) rather than blast tab output files. No bestHit approach will be taken into consideration.')
parser.add_argument('-d', '--hmmerFile', action='store_true', help='Process hmmscan output. The files should be on the parseable table output format (tblout or domtblout or pfamtblout arguments on hmmscan)')
parser.add_argument('-e', '--eggNOG', action='store_true', help='Set this argument if annotations were done against new eggNOG 4.1 version (implies -d).')
parser.add_argument('-b', '--bestHitOnly', action='store_true', help='Set this argument if you only want to process the best hit per gene.')
parser.add_argument('-a', metavar='addDescription', help='Two-column tab separated file giving the annotation description for each ID. ID must come in 1st col, descritption in 2nd.')
args = parser.parse_args()



if __name__ == "__main__":
###############################################################################
########################## COMMAND LINE PARSING ###############################
##############################################################################

    if len(args.input_files) < 2 :
        sys.stderr.write('Error: %s\n' % 'You should specify at least two annotation files.')
        sys.exit(2)

    if args.GIlist and args.hmmerFile:
        sys.stderr.write('Error: %s\n' % 'Please decide if input files are a ID list or a hmmer output file, not both.')
        sys.exit(2)

    elif args.GIlist and args.bestHitOnly:
        sys.stderr.write('Error: %s\n' % 'When gene identifiers files are provided (-g argument), no hits are being analysed. Please remove -b argument.')
        sys.exit(2)

    elif args.GIlist:
        dict_unique,dict_repeated=processFromIDs(args.input_files)
        intersection(dict_unique,dict_repeated,args.o, args.a)

    elif args.hmmerFile:
        annotated_genes,dict_unique,dict_repeated = processFromHmmer(args.input_files,args.bestHitOnly,args.eggNOG)
        intersection(dict_unique,dict_repeated,args.o,args.a)

    elif args.eggNOG:
        sys.stderr.write('Error: %s\n' % 'eggNOG v4.1 annotations are only available through an hmmscan search. Please set -d argument.')
        sys.exit(2)

    else:
        annotated_genes,dict_unique,dict_repeated = processFromBlastTab(args.input_files,args.bestHitOnly)
        intersection(dict_unique,dict_repeated,args.a)