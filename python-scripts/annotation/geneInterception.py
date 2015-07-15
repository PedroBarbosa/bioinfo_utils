__author__ = 'pedro'

import argparse
import sys
import os
import csv
import collections
def processFromBlastTab(inputFiles,bestHit):
    mydict_unique = {} #dictionary of all unique IDs obtained annotation file
    mydict_repeated = {} #dictionary of the repeated IDs per annotation file
    for filename in inputFiles:
        list_name = filename
        list_name_repeated = filename + "_repeated"
        mydict_unique[list_name] = []
        mydict_repeated[list_name_repeated] = []

        with open(filename) as file:
            previous_query = ""
            query = ""
            annotated_reads = 0
            print("Processing file " + filename + "..")
            if bestHit:

                for line in file:
                    if not line.startswith('#') and line.rstrip():

                        query = line.split()[1]
                        if query == previous_query:
                            previous_query = query

                        else:
                            annotated_reads += 1
                            hit = line.split("\t")[2]
                            values = mydict_unique[list_name]
                            if hit not in values:
                                values.append(hit)
                                mydict_unique[list_name] = values
                            else:
                                values_repeated = mydict_repeated[list_name_repeated]
                                values_repeated.append(hit)
                                mydict_repeated[list_name_repeated] = values_repeated




            else:
                for line in file:
                    if not line.startswith('#'):

                        query = line.split()[1]
                        if query != previous_query:
                            annotated_reads += 1

                        hit = line.split("\t")[2]
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
    mydict_unique = {} #dictionary of all unique IDs obtained annotation file
    mydict_repeated = {} #dictionary of the repeated IDs per annotation file
    for filename in inputFiles:


        list_name = filename
        list_name_repeated = filename + "_repeated"
        mydict_unique[list_name] = []
        mydict_repeated[list_name_repeated] = []

        with open(filename) as file:
            previous_query = ""
            query = ""
            annotated_reads = 0
            print("Processing file " + filename + "..")
            if bestHit:

                for line in file:
                    if not line.startswith('#') and line.rstrip():

                        query = line.split()[2]
                        if query == previous_query:
                            previous_query = query

                        else:
                            annotated_reads += 1

                            if eggNOG:
                                hit = line.split("\t")[0].split(".")[1] #eggnog ortholog group in an hmmscan hits file

                            else:

                                hit = line.split()[0]

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
                    if not line.startswith('#'):
                            query = line.split()[2]
                            if query != previous_query:
                                annotated_reads += 1

                            if eggNOG:
                                hit = line.split("\t")[0].split(".")[1] #eggnog ortholog group in an hmmscan hits file
                            else:

                                hit = line.split()[0]

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


def processFromIDs(inputFiles,bestHit):
    #createLists
    mydict_unique = {}
    mydict_repeated = {}
    list_all = []
    for filename in inputFiles:

        with open(filename) as file:
            list_name = filename
            list_name_repeated = filename + "_repeated"
            mydict_unique[list_name] = []
            mydict_repeated[list_name_repeated] = []
            for GI in file:
                values = mydict_unique[list_name]
                if GI not in values:
                    values.append(GI)
                    mydict_unique[list_name] = values
                else:
                    #update retepated IDs dict
                    values_repeated = mydict_repeated[list_name_repeated]
                    values_repeated.append(GI)
                    mydict_repeated[list_name_repeated] = values_repeated

    return mydict_unique,mydict_repeated

def intersection(dict_uniq,dict_repeat):
    #Get interception list
    print("\nCalculating interception..")
    print("Generating stats..\n")
    interception = dict_uniq.itervalues().next()
    previous_k = dict_uniq.iteritems().next()[0].split("/")[-1]
    max_possible = 100000000000
    for k,v in dict_uniq.iteritems():
        interception = list(set(interception) & set(v))
        if len(v) < max_possible:
            max_possible = len(v)
        #print("Annotations shared between %s and %s:\t%s" %  (previous_k, k.split("/")[-1],len(paired_interception)))
        #print(k.split("/")[-1] + "\t" + str(len(interception)))
    print ("Maximum number of possible overlaps (Size of the smallest set of IDs in all files):\t%s" % max_possible)
    print ("Number of annotations shared:\t%s (%s)" % (len(interception),round(float(len(interception)) / float(max_possible) * 100,2)))
    #Percentage
    for k,v in dict_uniq.iteritems():
        percentage = round(float(len(interception)) / float((len(v))) * 100,2)
        print("Percentage of %s shared across all annotations:\t%s" % (k.split("/")[-1], percentage))


    
    ##Unique analysis
    newDict={}
    print("\n")
    print("#################Single copy genes analysis ################\n######################################################")
    for k,v in dict_uniq.iteritems():

        for id in v:

            if id in newDict.keys():
                newDict[id].add(k.split("/")[-1])
            else:
                newDict[id] = set()
                newDict[id].add(k.split("/")[-1])

    outputFile = "uniqueIDs_output.txt"
    print("Writing output file to %s directory..\n\n" % os.path.abspath(outputFile))
    with open(outputFile, "w") as csvfile:
        writer = csv.writer(csvfile,dialect=csv.excel_tab)
        writer.writerow(("#List of unique IDs to each annotation",''))
        for id in sorted(newDict,key=newDict.get, reverse=True):
            print(id,newDict[id])
            if len(newDict[id]) == 1: #check if only one annotation has this id
                writer.writerow((''.join(newDict[id]),id))
    csvfile.close()



    #Repeated genes analysis
    newDict_rep = {}
    print("\n")
    print("####################Repeated genes analysis ################\n################################################################")
    for k,v in dict_repeat.iteritems():
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
    for id in sorted(newDict_rep, key=newDict_rep.get, reverse=True):
        print (id,newDict_rep[id])




parser = argparse.ArgumentParser(description='Script to check the interception of the genes present in different genome annotations (default blastTAB format).')
parser.add_argument(dest='input_files', metavar='annotated_files', nargs='+', help='Annotation files to be analyzed (minimum 2).')
parser.add_argument('-g', '--GIlist', action='store_true', help='Process gene identifiers (one per line) rather than a blast tab output file.')
parser.add_argument('-d', '--hmmerFile', action='store_true', help='Process hmmscan output. The file should be on the parseable table output format (tblout or domtblout or pfamtblout arguments on hmmscan)')
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

    elif args.GIlist:
        dict_unique,dict_repeated=processFromIDs(args.input_files, args.bestHitOnly)
        intersection(dict_unique,dict_repeated)

    elif args.hmmerFile:

        annotated_genes,dict_unique,dict_repeated = processFromHmmer(args.input_files,args.bestHitOnly,args.eggNOG)
        intersection(dict_unique,dict_repeated)

    elif args.eggNOG:
        sys.stderr.write('Error: %s\n' % 'eggNOG v4.1 annotations are only available through an hmmscan search. Please set -d argument.')
        sys.exit(2)

    else:
        annotated_genes,dict_unique,dict_repeated = processFromBlastTab(args.input_files,args.bestHitOnly)
        intersection(dict_unique,dict_repeated)