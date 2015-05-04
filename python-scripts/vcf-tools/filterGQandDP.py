__author__ = 'pedro'

import argparse
import os
import csv
import numpy as np
from collections import Counter
from collections import defaultdict
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt



def processFile(vcf):
    headers=[]
    map_DP=defaultdict(list)
    map_GQ=defaultdict(list)
    with open(vcf[0]) as file:
        for line in file:
            if line[0] != '#':
                columns = line.rstrip().split('\t')
                i=0
                genotypes = columns[9:]

                for sample in genotypes:

                    if "./." in sample:
                        map_DP[headers[i]].append("noSNP")
                        map_GQ[headers[i]].append("noSNP")

                    else:

                        #get values
                        dp=int(sample.split(':')[2])
                        gq=int(sample.split(':')[3])

                        #add new value to the dict
                        map_DP[headers[i]].append(dp)
                        map_GQ[headers[i]].append(gq)

                    i+=1

            elif line.startswith("#CHROM"):
                headers = line.rstrip().split('\t')[9:]

    file.close()
    return  map_DP,map_GQ,headers


def write2file(dictDP, dictGQ, headers,outputfile):

    if os.path.exists(outputfile):
        os.remove(outputfile)
    with open(outputfile, "w") as csvfile:
        writer = csv.writer(csvfile,dialect=csv.excel_tab)
        #print(headers)
        writer.writerow(('Total number of potential SNPs in file:', len(dictDP[headers[0]])))
        writer.writerow('')

        #Coverage analysis
        for k,v in dictDP.iteritems():

            #get only integers
            list_integers = [e for e in v if isinstance(e, int)]
            data = Counter(list_integers)
            writer.writerow(('Number of potential SNPS sites where sample ' + k + ' does not have have an alternative allele against the reference', v.count('noSNP')))
            writer.writerow(('Mean Depth of Coverage value for sample ' + k + ":", np.mean(list_integers)))
            writer.writerow(('Maximum value of Depth of Coverage for sample ' + k + ":", max(list_integers)))
            writer.writerow(('Minimum value of Depth of Coverage for sample ' + k + ":", min(list_integers)))
            writer.writerow(('Most common value of Depth of Coverage for sample ' + k + ":", data.most_common(1)))
            writer.writerow('')

        #Genotype quality analysis
        for k,v in dictGQ.iteritems():
            #get only integers
            list_integers2 = [e for e in v if isinstance(e, int)]
            data2 = Counter(list_integers2)
            writer.writerow(('Mean Genotype Quality value for sample ' + k + ":", np.mean(list_integers2)))
            writer.writerow(('Maximum value of Genotype Quality for sample ' + k + ":", max(list_integers2)))
            writer.writerow(('Minimum value of Genotype Quality for sample ' + k + ":", min(list_integers2)))
            writer.writerow(('Most common value of Genotype Quality for sample ' + k + ":", data2.most_common(1)))
            writer.writerow('')

    csvfile.close()



def drawPlot(dictDP,dictGQ,headers,outputfile):


    for key,value in dictDP.iteritems():
        int_list=[e for e in value if isinstance(e, int)]

        plt.hist(int_list,alpha=0.5)
        plt.xlabel('Value of Coverage')
        plt.ylabel('Counts')
        plt.title(r'Histogram of Depth of Coverage for the ' + key + " sample")


        # add a 'best fit' line
        #y = mlab.normpdf(bins, mu, sigma)
        #plt.plot(bins, y, 'r--')


        # Tweak spacing to prevent clipping of ylabel
        #plt.subplots_adjust(left=0.15)
        #plt.show()

parser = argparse.ArgumentParser(description='Script to analyze the genotype quality and depth of coverage values in a vcf file.')
parser.add_argument(dest='vcf_file', metavar='vcf_file', nargs=1,help='File to be processed.')
parser.add_argument(dest='output_file', metavar='output_file', nargs=1, help='Output file')
args = parser.parse_args()


if __name__ == "__main__":
    print("Processing vcf file..")
    dictDP, dictGQ, headers= processFile(args.vcf_file)
    print("Writing the output to " + args.output_file[0] + "..")
    write2file(dictDP, dictGQ, headers, args.output_file[0])
    print("Drawing plots..")
    drawPlot(dictDP,dictGQ,headers,args.output_file[0])

