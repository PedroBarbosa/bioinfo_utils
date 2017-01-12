import argparse


from lib2to3.pgen2.tokenize import double3prog
from decimal import *
from collections import OrderedDict
__author__ = 'pedro'

parser = argparse.ArgumentParser(description='This is a script to correct the evalue of rapsearch results')
parser.add_argument('-i', dest='blast_results_file', help='Blast results file', required=True)
parser.add_argument('-e', dest='evalue_threshold', help='Threshold to accept hits', required=True)
parser.add_argument('-o', dest='output_file', help='Corrected evalue file', required=True)
args = parser.parse_args()


def processEvalue(blastFile, eval_threshold):
    dict_final=OrderedDict()
    features_with_hits=0
    previous_feature=""
    print("Reading and processing file ...")
    with open(args.output_file,'w') as outfile:
        with open(blastFile) as file:
            for line in file:
                if not line.startswith('#'):
                    vector = line.split("\t") #read itself
                    try:
                        #eValue = round(10 ** (Decimal(vector[10])),8) #new eValue
                        eValue = pow(10,float(vector[10])) 
                        if previous_feature != vector[0]:
                            features_with_hits+=1
                            #print(features_with_hits)
                        previous_feature=vector[0]
                        
                        if eValue <= float(eval_threshold):
                            vector[10] = format(eValue,'2e')
                            if vector[0] in dict_final:
                                dict_final[vector[0]].append(tuple(vector[1:]))
                            else:
                                dict_final[vector[0]] = [tuple(vector[1:])]
                    except:
                        pass   
                else:
                    vector = line.split("\t") #read itself
                    if len(vector) > 10:
                        vector[10] = "e-value"
                        outfile.write('\t'.join(vector))
                    else:
                        outfile.write(line)
                        
        for feature, hits in dict_final.items():
            for tup in hits:
                outfile.write(feature + "\t" + "\t".join(x for x in tup))

    print("Number of features with hits in the original file\t%i" % features_with_hits)
    print("Number of features kept after processing e-values\t%i" % len(dict_final))
    print("Done.")

processEvalue(args.blast_results_file, args.evalue_threshold)
