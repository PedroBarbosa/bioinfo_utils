__author__ = 'pedro'

import argparse
import logging
from collections import OrderedDict
import os

def createDict(fileList):
    dict= OrderedDict()
    listStats=[]

    key=""
    for file in fileList:
        strategy = os.path.basename(file).rsplit('.',1)[0]

        with open(file, 'r') as infile:
            previous=False
            for line in infile:
                if line.startswith("Generating") and previous:

                    if key in dict:
                        dict[key].append(listStats)
                    else:
                        dict[key] = []
                        dict[key].append(listStats)

                    listStats=[]
                    key=line.split(" ")[4].split("/")[-1]
                    listStats.append(strategy)

                elif line.startswith("Generating") and not previous: #if first line
                    key=line.split(" ")[4].split("/")[-1]
                    listStats.append(strategy)
                    previous=True

                elif "Number" in line:
                    number =[int(s) for s in line.split() if s.isdigit()]
                    listStats.append(number[0])

            #last mapping stat
            if key in dict:
               dict[key].append(listStats)
            else:
                dict[key] = []
                dict[key].append(listStats)
                listStats = []
        infile.close()

    return  dict

def writeOutput(dict):

    with open(os.getcwd() + "/mappingComparison.xlsx", "w") as outfile:



        dict_stat = OrderedDict()
        i=1
        for value in dict.values():
            for strategy in value:
                if strategy[0] not in dict_stat:
                    dict_stat[strategy[0]] = i
                    i += 1

        aux_list = [''] * (len(dict_stat) - 1)

        outfile.write(""+ "\t" + "#Alignments" + "\t" + '\t'.join(aux_list) + "\t" + "#Unmapped reads" + '\t'+ '\t'.join(aux_list) + "\t" + "#Primary and linear alignments"
                      + '\t'+ '\t'.join(aux_list) + "\t" + "Primary and linear mapq > 10" + '\t'+ '\t'.join(aux_list) + "\t" + "#Secondary alignments" + '\t'+'\t'.join(aux_list)
                      + "\t" + "#Chimeric alignments" + '\t'+'\t'.join(aux_list) + "\t" + "#Proper pair" + '\t'+'\t'.join(aux_list) + "\t" +"#Proper reads FR" + "\t" +'\t'.join(aux_list) + '\t' +
                      "#Proper reads RF" +"\t" +'\t'.join(aux_list) + "\t" + "#Unique" +'\t'+ '\t'.join(aux_list)+ '\n')

        header=False
        for key,value in dict.items():
            strategies,alignments,unmapped,primary_linear,primary_q10,secondary,chimeric,proper,proper_FR,proper_RF,unique = [],[],[],[],[],[],[],[],[],[],[]

            print(key,value)
            for strategy in value:
                strategies.append(str(strategy[0]))
                alignments.append(str(strategy[1]))
                unmapped.append(str(strategy[2]))
                primary_linear.append(str(strategy[3]))
                primary_q10.append(str(strategy[4]))
                secondary.append(str(strategy[5]))
                chimeric.append(str(strategy[6]))
                proper.append(str(strategy[7]))
                proper_FR.append(str(strategy[8]))
                proper_RF.append(str(strategy[9]))
                unique.append(str(strategy[10]))


            if header:
                outfile.write(key + '\t' + '\t'.join(alignments) + '\t' + '\t'.join(unmapped) + '\t' + '\t'.join(primary_linear) + '\t' + '\t'.join(primary_q10) + '\t'
                              + '\t'.join(secondary) + '\t' + '\t'.join(chimeric) + '\t' + '\t'.join(proper) + '\t' + '\t'.join(proper_FR) + '\t'
                              + '\t'.join(proper_RF) + '\t' + '\t'.join(unique) + '\t' +'\n')
            else:
                outfile.write("Sample" '\t' + '\t'.join(strategies) + '\t' + '\t'.join(strategies) + '\t' + '\t'.join(strategies) + '\t' + '\t'.join(strategies) + '\t'
                              + '\t'.join(strategies) + '\t' + '\t'.join(strategies) + '\t' + '\t'.join(strategies) + '\t' + '\t'.join(strategies) + '\t'
                              + '\t'.join(strategies) + '\t' + '\t'.join(strategies) + '\t' +'\n')
                header=True
                outfile.write(key + '\t'+'\t'.join(alignments) + '\t' + '\t'.join(unmapped) + '\t' + '\t'.join(primary_linear) + '\t' + '\t'.join(primary_q10) + '\t'
                              + '\t'.join(secondary) + '\t' + '\t'.join(chimeric) + '\t' + '\t'.join(proper) + '\t' + '\t'.join(proper_FR) + '\t'
                              + '\t'.join(proper_RF) + '\t' + '\t'.join(unique) + '\t' +'\n')
      #      if len(value) == 2:# == len(files):
      #          logging.error("%s does not have statistics in all (%s) files provided. Can't generate excel." % (key,len(files)))
      #      else:
      #          continue
def main():
    parser = argparse.ArgumentParser(description='Script to generate excel tables which compare different outputs of "generatingMappingStats.sh" script. At least 2 files are required.')
    parser.add_argument(dest='txt_files', metavar='statsFiles', nargs='+', help='Stats files to process. [minimum 2].')
    parser.add_argument('-i', '--ignore', action="store_true", help='Write excel file, if only one input file provided.')
    args = parser.parse_args()

    if len(args.txt_files) < 2 and not args.ignore:
        logging.error("Please provide at least 2 files to compare.")
        exit(1)
    elif len(args.txt_files) >= 2 and args.ignore:
        logging.error("You provided more than one input file. No need to set '-i' argument.Exiting..")
        exit(1)
    elif args.ignore:
        logging.info("One file provided. Excel will be generated for this file.")
        dict = createDict(args.txt_files)
        writeOutput(dict)
    elif len(args.txt_files) >= 2:
        logging.info("Correct input. Generating excel.")
        dict = createDict(args.txt_files)
        writeOutput(dict)
    else:
        logging.error("Some error happened.")
        exit(1)

if __name__ == "__main__":
    main()