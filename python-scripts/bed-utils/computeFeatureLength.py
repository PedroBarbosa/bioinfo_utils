import argparse
import pybedtools
from pybedtools import BedTool
from pybedtools import featurefuncs
import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO,  format='%(asctime)s %(message)s')
import matplotlib.pyplot as plt
import seaborn as sns

def addIntronLengthToScoreCol(feature):
    feature[4] = len(feature)
    if "_" in feature[3]:
        feature[3] = feature[3].split("_")[1]
    return feature

def filterByLength(feature,length):
    if len(feature) < length:
        return feature

def plotIntrons(bedobj,outbasename,featureType,length=0):
    df = bedobj.to_dataframe(names=['chrom', 'start', 'stop', 'name', featureType + '_length', 'strand'])
    plt.figure(figsize=(10, 9))
    sns.stripplot(data=df,y=featureType + "_length", x="name",jitter=True)
    plt.xticks(rotation='vertical',fontsize=12)
    plt.yticks(fontsize=14)
    plt.xlabel("")
    plt.ylabel(featureType + " length (bp)", fontsize=18)
    if length == 0:
        plt.ylim(0, df[featureType + '_length'].max() + 1000)
        #plt.title("Full intron")
        plt.savefig(outbasename + ".pdf", format='pdf')
    else:
        plt.ylim(0, length)
        plt.title("{} bp {} size threshold".format(length,featureType),fontsize=20)
        plt.savefig(outbasename + "_" + str(length) + ".pdf", format='pdf')

def processIntronicBed(inbed,lengths,outbasename,featureType):
    logging.info("Creating BedTool object.")
    myBedTool = BedTool(inbed)
    logging.info("Writing " + featureType + " length.")
    bedobj=myBedTool.each(addIntronLengthToScoreCol).saveas(outbasename + ".bed")
    plotIntrons(bedobj,outbasename,featureType,0)
    for l in lengths:
        smallIntObj=bedobj.each(filterByLength,l).saveas(outbasename + "_smaller" + str(l) + ".bed")
        plotIntrons(smallIntObj,outbasename,featureType,l)

def main():
    parser = argparse.ArgumentParser(description='Script to compute size of the features as well as plot their distribution based on the bed feature name')
    parser.add_argument(dest='bed', help='Path to the bed file')
    parser.add_argument('-l','--feature_length_bins', nargs='+', required=True, type=int,help='Intron size thresholds to write down and plot smaller introns')
    parser.add_argument("-o", "--output", required=True, help='Basename to write output files.')
    parser.add_argument("-t","--featureType",required=True, choices=("exon","intron"),help='Features represented in the input file')
    args = parser.parse_args()

    processIntronicBed(args.bed, args.feature_length_bins,args.output,args.featureType)

if __name__ == "__main__":
    main()
