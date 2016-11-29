import argparse
import logging
import sys
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
import urllib.request
import json
import collections
from collections import defaultdict
import os.path

def nogs2list(nogsFile, tab, col):
    listnogs=[]
    with open(nogsFile, 'r') as file:
        if tab:
            for line in file:
                if not line.startswith("#"):
                    nog_line = line.rstrip()
                    if len(nog_line.split("\t")) == 1:
                        logging.error("Your input file in single column. Please remove the '-tab' and '-col' arguments from your command.")
                        exit(1)

                    elif len(nog_line.split("\t")) < col:
                        print(len(nog_line.split("\t")))
                        logging.error("Tab separated file doesn't have column number you supplied in '-col' argument. Please change the column number argument.")
                        exit(1)
                    nog = nog_line.split("\t")[col-1]

                    if not nog in listnogs and not "available" in nog:
                        listnogs.append(nog.rstrip())
                    elif not "available" in nog:
                        logging.warning("%s id repeated in file" % nog)
        else:
            for line in file:
                if not line.startswith("#"):
                    nog = line.rstrip()
                    if not nog in listnogs:
                        listnogs.append(nog.rstrip())
                    else:
                        logging.warning("%s id repeated in file" % nog)

    return listnogs


def getAnnotationFromAPI(listNOGs):
    logging.info("Using EggNOG API to retrive NOGs annotations..")
    #http://eggnogapi.embl.de/nog_data/[dataformat]/[attributes]/[nogname]

    #[dataformat] currently supported data formats are:
    #    json: data is returned in compressed JSON format (allows for multiple attributes)
    #    text: data is returned as a text string
    #    file: data is returnend as text file or as tar.gz (allows for multiple attributes)
    #[attributes] One or more (comma separated):
    #    fasta, raw_alg, trimmed_alg, tree, go_terms, domains
    #[nogname] A valid EggNOG group name (i.e. ENOG410ZSWV or COG0575)

    lessFreqGO = 0
    lessFreqDom = 0
    out_dict = collections.OrderedDict()

    f = urllib.request.urlopen('http://eggnogapi.embl.de/nog_data/json/domains,go_terms/ENOG41010CW').read()
    data = json.loads(f.decode('utf-8'))

    i=1
    ##Assumes that domain dict comes first in the iteration

    for nog in listNOGs:
        logging.info("Processing %s NOG.. %i" % (nog,i))
        i+=1
        f = urllib.request.urlopen('http://eggnogapi.embl.de/nog_data/json/domains,go_terms/'+nog).read()
        data = json.loads(f.decode('utf-8'))

        for k,subdict in iter(data.items()):

            if k == "domains":
                if not nog in out_dict.keys():
                    out_dict[nog] = ["-"]*5

                index = 0
                for domain_database,listDomainsPerDb in iter(subdict.items()):
                    if domain_database == "SMART":
                        index = 0
                    elif domain_database == "PFAM":
                        index = 1
                    domain_value = ""
                    oneHighFreqDomain = False
                    for domain in listDomainsPerDb:
                        if float(domain[2]) > 90: #if frequency of this domain in the sequences used to produced this NOG is higher than 90%
                            domain_value = domain_value + domain[0] + "; "
                            oneHighFreqDomain = True
                        else:
                            lessFreqDom += 1
                    if oneHighFreqDomain:
                        out_dict[nog][index] = domain_value
                    else:
                        out_dict[nog][index] = "No domain found with more than 90% frequency"



            elif k == "go_terms":
                if not nog in out_dict.keys():
                    out_dict[nog] = ["-"]*5

                index = 0
                for go_category, gosListPerCategory in iter(subdict.items()):
                    #print(len(subdict))
                    if go_category == 'Cellular Component':
                        #print("ola")
                        index = 2
                    elif go_category == 'Molecular Function':
                        #print("ola2")
                        index = 3
                    elif go_category == 'Biological Process':
                        #print("ola3")
                        index = 4

                    go_value = ""
                    for go in gosListPerCategory:
                        if float(go[4]) > 90: #if frequency of this GO in the sequences used to produced this NOG is higher than 90%
                            go_value = go_value + go[0] + " [" + go[1] + "]" + "; "
                        else:
                            lessFreqGO += 1
                    if go_value == "":
                        out_dict[nog][index] = "No GO found with more than 90% frequency."
                    else:
                        out_dict[nog][index] = go_value

    logging.info("Number of NOGs with retrieved annotation\t%i" % len(out_dict))
    return out_dict


def addCOGcategoriesFromFile(aux_file,out_dict):
    logging.info("Adding COG categories and descriptions to output..")
    nogs_with_match = 0
    with open(aux_file, 'r') as cogFile:
        for line in cogFile:
            nog = line.rstrip().split("\t")[1]
            if line.rstrip().split()[1] in out_dict.keys():
                nogs_with_match += 1
                cog_category = line.rstrip().split("\t")[4]
                description = line.rstrip().split("\t")[5]
                if description == "NA":
                    out_dict[nog].extend([cog_category, "-"])
                else:
                    out_dict[nog].extend([cog_category, description])
    cogFile.close()
    percentage = round(float(nogs_with_match/len(out_dict) * 100),4)
    logging.info("Number of NOGs for which additional information was added from the COG categories file\t%i(%s%%)", nogs_with_match,str(percentage))
    return out_dict

def writeOutput(outFile, out_dict, addAnn):
    logging.info("Writing to output file..")
    with open(outFile, 'w') as out_file:
        if addAnn:
            out_file.write("#NOG" + "\t" + "SMART domains" + "\t" + "PFAM_domains" + "\t" + "GOs Cellular Componet" + "\t" + "GOs Molecular Function" + "\t" +
        "GOs Biological Process" + "\t" + "COG category" + "\t" + "COG description" + "\n")
        else:
            out_file.write("#NOG" + "\t" + "SMART domains" + "\t" + "PFAM_domains" + "\t" + "GOs Cellular Componet" + "\t" + "GOs Molecular Function" + "\t" +
        "GOs Biological Process" + "\n")
        for nog, annotations in iter(out_dict.items()):
            #print(nog,annotations)
            out_file.write(nog + "\t" + "\t".join(i for i in annotations) + "\n")


def writeOutputFromTabSeparatedFile(outFile, out_dict, addAnn, inputFile,col):
    with open(inputFile, 'r') as file:
        with open(outFile, 'w') as out_file:
            for line in file:

                if not line.startswith("#"):
                    numb_col = len(line.split())
                    aux_list = ['']*numb_col
                    break
            file.close()
            if addAnn:
                out_file.write("#" + "\t".join(aux_list) + "SMART domains" + "\t" + "PFAM_domains" + "\t" + "GOs Cellular Componet" + "\t" + "GOs Molecular Function" + "\t" +
        "GOs Biological Process" + "\t" + "COG category" + "\t" + "COG description" + "\n")
            else:
                out_file.write("#" + "\t".join(aux_list) + "SMART domains" + "\t" + "PFAM_domains" + "\t" + "GOs Cellular Componet" + "\t" + "GOs Molecular Function" + "\t" +
        "GOs Biological Process" + "\n")
            with open(inputFile,'r') as infile:
                for line in infile:
                    if not line.startswith("#"):
                        nog = line.split()[col-1]
                        if nog in out_dict.keys():
                            out_file.write(line.rstrip() + '\t' + '\t'.join(out_dict[nog]) + "\n")
                        else:
                            out_file.write(line.rstrip() + "\n")

def writeOutputFromHmmer(outFile, out_dict, addAnn, nog2genedict):
    logging.info("Writing to output file..")
    with open(outFile, 'w') as out_file:
        if addAnn:
            out_file.write("#Feature" + "\t" + "NOG" + "\t" + "SMART domains" + "\t" + "PFAM_domains" + "\t" + "GOs Cellular Componet" + "\t" + "GOs Molecular Function" + "\t" +
        "GOs Biological Process" + "\t" + "COG category" + "\t" + "COG description" + "\n")
        else:
            out_file.write("#Feature" + "\t" + "NOG" + "\t" + "SMART domains" + "\t" + "PFAM_domains" + "\t" + "GOs Cellular Componet" + "\t" + "GOs Molecular Function" + "\t" +
        "GOs Biological Process" + "\n")
        for nog, annotations in iter(out_dict.items()):
            #print(nog,annotations)
            if nog in nog2genedict:
                if len(nog2genedict[nog]) > 1 :
                    logging.info("Multiple genes are assigned to %s NOG" % nog)
                    val = nog2genedict[nog]
                    for v in val:
                        out_file.write(v + "\t" + nog + "\t" + "\t".join(i for i in annotations) + "\n")
                 else:
                      out_file.write(str(nog2genedict[nog]) + "\t" + nog + "\t" + "\t".join(i for i in annotations) + "\n")
            else:
                logging.error("Some weird error happened. Contact Pedro.")
                exit(1)


def processFromHmmer(nogs):

    listNogs=set()
    dict=defaultdict(list)
    with open(nogs,'r') as infile:
        annotated_genes = 0
        for line in infile:

            previous_query = ""
            print("Processing hmmer parsable file " + nogs + "..")

            if not line.startswith('#') and line.rstrip():
                query = line.split()[2]
                if query != previous_query:
                    annotated_genes += 1
                    hit = line.split("\t")[0].split(".")[1] #eggnog ortholog group in an hmmscan hits file
                    previous_query = query
                    listNogs.add(hit)
                    dict[hit].append(query)

        print("%s annotated genes." % (annotated_genes))
    return listNogs,dict

def main():
    parser = argparse.ArgumentParser(description='Script add annotation information of an eggNOG record. Uses eggNOG API to retrieve data. Outputs a tab separarted file.')
    parser.add_argument(dest='nogs', metavar='nogsFile', help='File listing the NOGs to process. One perl line.')
    parser.add_argument(dest='output', metavar='outputFile', help='Output File will be a tab delimited file to open in excel or related programs.')
    parser.add_argument('-add' ,metavar ='-addNOGAnnotationFile', help='Additional file with information about each NOG, particularly COG categories and descriptions. Downloadable from eggNOG website.')
    parser.add_argument('-hmmer', metavar = '--hmmerFile',action='store_true', help="NOGs file is a hmmer table output file. If set '-tab' and '-col' should not be set.")
    parser.add_argument('-tab', '--tabSeparated', action='store_true', help="Flag indicating that NOGs file is tab separated. If set, requires '-col' to be set too.")
    parser.add_argument('-col', metavar = '-column', type=int, nargs=1, help="Column number in the tab separated file in which NOG records are displayed.")
    args = parser.parse_args()

    if args.tabSeparated and not args.col:
        logging.error("If -tab is provided, please set also '-col' argument.")
        exit(1)
    elif args.col and not args.tabSeparated:
        logging.error("Please set '-col' argument only if input file in tab separated. Either remove '-col' or add '-tab' to the command.")
        exit(1)

    if args.hmmer and args.tabSeparated or args.hmmer and args.col:
        logging.error("-tab or -col arguments not available when -hmmer flag is set")
        exit(1)


    if args.hmmer:
        listNOGs, nog2geneDict = processFromHmmer(args.nogs)
        api_dict = getAnnotationFromAPI(listNOGs)
        if args.add:
            if os.path.isfile(args.add):
                out_dict = addCOGcategoriesFromFile(args.add,api_dict)
            else:
                logging.error("Please provide a valid file in the '-add' argument.")
                exit(1)
        writeOutputFromHmmer(args.output,api_dict,args.add,nog2geneDict)


    else:
        listNOGs = nogs2list(args.nogs, args.tabSeparated, args.col[0])
        out_dict = getAnnotationFromAPI(listNOGs)
        if args.add:
            if os.path.isfile(args.add):
                out_dict = addCOGcategoriesFromFile(args.add,out_dict)
            else:
                logging.error("Please provide a valid file in the '-add' argument.")
                exit(1)
        if args.tabSeparated:
            writeOutputFromTabSeparatedFile(args.output,out_dict, args.add, args.nogs,args.col[0])
        else:
            writeOutput(args.output,out_dict, args.add)



if __name__ == "__main__":
    main()