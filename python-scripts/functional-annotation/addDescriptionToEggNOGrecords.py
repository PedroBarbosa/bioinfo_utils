import argparse
import logging
import sys
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
import urllib.request
import json
from collections import defaultdict
import os.path

def nogs2list(nogsFile):
    listnogs=[]
    with open(nogsFile, 'r') as file:
        for nog in file:
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
    out_dict = defaultdict()

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





def main():
    parser = argparse.ArgumentParser(description='Script add annotation information of an eggNOG record. Uses eggNOG API to retrieve data. Outputs a tab separarted file.')
    parser.add_argument(dest='nogs', metavar='nogsFile', help='File listing the NOGs to process. One perl line.')
    parser.add_argument(dest='output', metavar='outputFile', help='Output File will be a tab delimited file to open in excel or related programs.')
    parser.add_argument('-add',metavar ='-addNOGAnnotationFile', help='Additional file with information about each NOG, particularly COG categories and descriptions. Downloadable from eggNOG website.')
    args = parser.parse_args()


    listNOGs = nogs2list(args.nogs)
    out_dict = getAnnotationFromAPI(listNOGs)
    if args.add:
        if os.path.isfile(args.add):
            out_dict = addCOGcategoriesFromFile(args.add,out_dict)
        else:
            logging.error("Please provide a valid file in the '-add' argument.")
            exit(1)
    writeOutput(args.output,out_dict, args.add)

if __name__ == "__main__":
    main()