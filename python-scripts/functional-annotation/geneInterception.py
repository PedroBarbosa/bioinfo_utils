__author__ = 'pedro'

import argparse
import sys
import os
import csv
from itertools import combinations
from collections import defaultdict
import subprocess

def processFromBlastTab(inputFiles,bestHit):
    mydict_unique = {} #dictionary of all unique IDs obtained annotation file
    mydict_repeated = {} #dictionary of the repeated IDs per annotation file
    for filename in inputFiles:
        list_name = filename
        list_name_repeated = filename + "_repeated"
        mydict_unique[list_name] = []
        mydict_repeated[list_name_repeated] = []
        annotated_reads = 0

        with open(filename) as file:

            previous_query = ""
            print("Processing file " + filename + "..")
            if bestHit:

                for line in file:
                    if len(line.split()) == 1:
                        print("%s file represents gene identifiers, please set '-g' argument accordingly. Exiting.." % filename)
                        exit(2)
                    if not line.startswith('#') and line.rstrip():

                        query = line.split()[0]
                        if "WP_" in query:
                            query = query.split("|")[-1]
                            print (query)
                        if query == previous_query:
                            previous_query = query

                        else:
                            annotated_reads += 1
                            hit = line.split()[1]
                            values = mydict_unique[list_name]
                            if hit not in values:
                                values.append(hit)
                                mydict_unique[list_name] = values
                            else:
                                values_repeated = mydict_repeated[list_name_repeated]
                                values_repeated.append(hit)
                                mydict_repeated[list_name_repeated] = values_repeated

                        previous_query = query


            else:
                for line in file:
                    if len(line.split()) == 1:
                        print("%s file represents gene identifiers, please set '-g' argument accordingly. Exiting.." % filename)
                        exit(2)

                    if not line.startswith('#'):

                        query = line.split()[0]
                        if query != previous_query:
                            annotated_reads += 1

                        hit = line.split()[1]
                        values = mydict_unique[list_name]
                        if hit not in values:
                            values.append(hit)
                            mydict_unique[list_name] = values
                        else:
                            values_repeated = mydict_repeated[list_name_repeated]
                            values_repeated.append(hit)
                            mydict_repeated[list_name_repeated] = values_repeated

                        previous_query = query

        print("%s annotated genes/reads in file %s!" % (annotated_reads,filename))
    return annotated_reads,mydict_unique,mydict_repeated


def processFromHmmer(inputFiles,bestHit,eggNOG):
    mydict_unique = defaultdict(list) #dictionary of all unique IDs obtained annotation file
    mydict_repeated = defaultdict(list) #dictionary of the repeated IDs per annotation file
    for filename in inputFiles:
        list_name = filename
        list_name_repeated = filename + "_repeated"
        mydict_unique[list_name] = []
        mydict_repeated[list_name_repeated] = []
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
                            mydict_unique[list_name].append(hit) if hit not in values else mydict_repeated[list_name_repeated].append(hit)
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
                        mydict_unique[list_name].append(hit) if hit not in values else mydict_repeated[list_name_repeated].append(hit)
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
            list_name_repeated = filename + "_repeated"
            mydict_unique[list_name] = []
            mydict_repeated[list_name_repeated] = []
            for gi in file:
                values = mydict_unique[list_name]
                GI = gi.rstrip()
                if GI not in values:
                    mydict_unique[list_name].append(GI)
                else:
                    mydict_repeated[list_name_repeated].append(GI)

    return mydict_unique,mydict_repeated

def intersection(dict_uniq,dict_repeat,outputBasename):
    #Get interception list
    print("\nCalculating interception..")
    print("Generating stats..\n")

    interception = set()

    for v in iter(dict_uniq.values()):
        interception.update(v)

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
        writer = csv.writer(csvfile,dialect=csv.excel_tab)
        i = 0
        listKeys = []
        dict_withIndex = {}
        for k in dict_uniq.keys():
            key = k.split("/")[-1]
            listKeys.append(key)
            dict_withIndex[key] = i
            i += 1


        writer.writerow(("#List of unique IDs to each annotation",'\t'.join(listKeys)))
        for k in sorted(newDict, key=lambda k: len(newDict[k]), reverse=True):
            final_list = ['no'] * len(listKeys)
            for annotations in newDict[k]:
                index = dict_withIndex[annotations]
                final_list[index] = 'yes'
            writer.writerow((k,'\t'.join(final_list)))

        #removeChar(outputBasename + "_interseptionTable.txt")

    csvfile.close()

    #Repeated genes analysis
    newDict_rep = {}
    print("\n")
    print("####################Repeated genes analysis ################\n################################################################")
    for k,v in iter(dict_repeat.items()):
        print("Number of genes with more than one copy in the genome for %s sample:\t%s" % (k.split("/")[-1].split("_")[:-1],len(set(v))))
        for repeated in v:

            if repeated in newDict_rep:
                if k.split("/")[-1] in newDict_rep[repeated].keys():
                    newDict_rep[repeated][k.split("/")[-1]] += 1
                else:
                    newDict_rep[repeated][k.split("/")[-1]] = 2

            else:
                newDict_rep[repeated] = {}
                newDict_rep[repeated][k.split("/")[-1]] = 2 #2nd copy because firts was not added to this dict previously (was added to the unique set)



    print("\nFrequency of the > 1x copy genes across samples:")
    if bool(newDict_rep):
        for id in sorted(newDict_rep, key=newDict_rep.get, reverse=True):
            print (id,newDict_rep[id])
    else:
        print("No genes found on this condition.")



#def removeChar(filename):
#    print(filename)
#    subprocess.call(['sed', '-i', 's/\"//g', filename])

parser = argparse.ArgumentParser(description='Script to check the interception of the genes present in different genome annotations (default blastTAB format).')
parser.add_argument(dest='input_files', metavar='annotated_files', nargs='+', help='Annotation files to be analyzed (minimum 2).')
parser.add_argument('-o', metavar='outputBasename', required = True, help='Output basename for which the output files will be written.')
parser.add_argument('-g', '--GIlist', action='store_true', help='Process gene identifiers (one per line) rather than blast tab output files. No bestHit approach will be taken into consideration.')
parser.add_argument('-d', '--hmmerFile', action='store_true', help='Process hmmscan output. The files should be on the parseable table output format (tblout or domtblout or pfamtblout arguments on hmmscan)')
parser.add_argument('-e', '--eggNOG', action='store_true', help='Set this argument if annotations were done against new eggNOG 4.1 version (implies -d).')
parser.add_argument('-b', '--bestHitOnly', action='store_true', help='Set this argument if you only want to process the best hit per gene.')
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
        intersection(dict_unique,dict_repeated,args.o)

    elif args.hmmerFile:
        annotated_genes,dict_unique,dict_repeated = processFromHmmer(args.input_files,args.bestHitOnly,args.eggNOG)
        intersection(dict_unique,dict_repeated,args.o)

    elif args.eggNOG:
        sys.stderr.write('Error: %s\n' % 'eggNOG v4.1 annotations are only available through an hmmscan search. Please set -d argument.')
        sys.exit(2)

    else:
        annotated_genes,dict_unique,dict_repeated = processFromBlastTab(args.input_files,args.bestHitOnly)
        intersection(dict_unique,dict_repeated)