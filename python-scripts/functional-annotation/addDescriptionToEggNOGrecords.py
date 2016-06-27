import argparse
import logging
import sys
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
import urllib.request
import json
from collections import defaultdict


def nogs2list(nogsFile):
    listnogs=[]
    with open(nogsFile, 'r') as file:
        for nog in file:
            if not nog in listnogs:
                listnogs.append(nog.rstrip())
            else:
                logging.warning("%s id repeated in file" % nog)
    return listnogs


def getAnnotationFromAPI(listNOGs,output):
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

    f = urllib.request.urlopen('http://eggnogapi.embl.de/nog_data/json/domains,go_terms/ENOG4100E81').read()
    data = json.loads(f.decode('utf-8'))
    print(data)
    i=1
    for nog in listNOGs:
        logging.info("Processing %s NOG.. %i" % (nog,i))
        i+=1
        f = urllib.request.urlopen('http://eggnogapi.embl.de/nog_data/json/domains,go_terms/'+nog).read()
        data = json.loads(f.decode('utf-8'))

        for k,subdict in iter(data.items()):
            out_dict[nog] = []
            if k == "go_terms":
                index = 0
                print(subdict)
                for go_category, gosListPerCategory in iter(subdict.items()):

                    if go_category == 'Molecular Funcion':
                        index = 0
                    elif go_category == 'Biological Process':
                        index = 1
                    elif go_category == 'Cellular Component':
                        index = 2

                    go_value = ""
                    for go in gosListPerCategory:
                        if float(go[4]) > 90: #if frequency of this GO in the sequences used to produced this NOG is higher than 90%
                            go_value = go_value + go[0] + "[" + go[1] + "]" + ";"
                        else:
                            lessFreqGO += 1
                    if go_value == "":
                        out_dict[nog].append("-")
                    else:
                        out_dict[nog].append(go_value)

                for k,v in iter(out_dict.items()):
                    print(k,v)


                 #print(k,subdict)
                 #print("\n\n\n")
    #         if k == "domains":
    #             for subkey,value in iter(data.items()):
    #                 if subkey != "dom_header":
    #
    #                     oneHighFreqDomain = False
    #                     out_dict[nog] = []
    #                     for domain_database,listDomainsPerDb in iter(value.items()):
    #                         domain_value = domain_database
    #                         for domain in listDomainsPerDb:
    #                             if float(domain[2]) > 90: #if frequency of this domain in the sequences used to produced this NOG is higher than 90%
    #                                 domain_value = domain_value + "**" + domain[0]
    #                                 oneHighFreqDomain = True
    #                             else:
    #                                 lessFreqDom += 1
    #                         if oneHighFreqDomain:
    #                             out_dict[nog].append(domain_value)
    #
    #                     for nog,domains in iter(out_dict.items()):
    #                         if len(domains) == 0:
    #                             out_dict[nog].extend(["No SMART domains", "No PFAM domains"])
    #                         elif len(domains) == 1 and any("PFAM" in s for s in domains):
    #                             domains.insert(0,"No SMART domains")
    #                             out_dict[nog] = domains
    #                         elif len(domains) == 1 and any("SMART" in s for s in domains):
    #                             out_dict[nog].append("No PFAM domains")
    #
    # for k,v in iter(out_dict.items()):
    #     print(k,v)


def main():
    parser = argparse.ArgumentParser(description='Script add annotation information of an eggNOG record. Uses eggNOG API to retrieve data. Outputs a tab separarted file.')
    parser.add_argument(dest='nogs', metavar='nogsFile', help='File listing the NOGs to process. One perl line.')
    parser.add_argument(dest='output', metavar='outputFile', help='Output File')
    parser.add_argument('-a',metavar ='-addNOGAnnotationFile', help='Additional file with information about each NOG, particularly COG categories and descriptions. Downloadable from eggNOG website.')
    args = parser.parse_args()


    listNOGs = nogs2list(args.nogs)
    getAnnotationFromAPI(listNOGs, args.output)

if __name__ == "__main__":
    main()