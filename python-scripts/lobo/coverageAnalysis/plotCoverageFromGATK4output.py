import argparse
import os
import sys
import logging
import glob
import subprocess
logging.basicConfig(stream=sys.stdout, level=logging.INFO,  format='%(asctime)s %(message)s')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from sklearn import preprocessing
from scipy.cluster import hierarchy
from scipy.spatial import distance
def wccount(filename):
    out = subprocess.Popen(['wc', '-l', filename],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT
                         ).communicate()[0]
    return int(out.partition(b' ')[0])

def countOcurrences(filename,char):
    with open(filename) as f:
        return sum(char in line for line in f)

#def plotWGShisto(filename):
def processTargetedExperiment(allMetrics,perTargetMetrics):
    logging.info("Starting targeted experiment analysis.")
    if len(list(glob.glob('*HS_metrics.txt'))) > 0 #if individual metrics files exist
        logging.info("We detected individual sample metrics files. Will check the low quality base pair distribution in the data")
        processIndividualTargetSamples()

    targetedAllMetricsProcess(allMetrics)

def processIndividualTargetSamples():
    dict_bq={}
    nleft=len(list(glob.glob('*HS_metrics.txt')))
    for indiv_metrics in list(glob.glob('*HS_metrics.txt')):
        sname=indiv_metrics.split("_")[0]
        logging.info("{} file! {} left..".format(indiv_metrics, nleft))
        nleft -= 1
        with open(indiv_metrics, "r") as f:
            baseQualities=[]
            for line in f:
                if "coverage_or_base_quality" in line:
                    for i, line in enumerate(f):
                        if i < 42:
                            baseQualities.append(int(line.rstrip().split("\t")[2]))
                        else:
                            if sname in dict_bq.keys():
                                logging.error("Repeated {} sample in directory.".format(indiv_metrics))
                                exit(1)
                            else:
                                dict_bq[sname] = baseQualities
                                break
        f.close()
    boxplot_dic={}
    for sample, baseQualities in dict_bq.items():
        totalbp=sum(baseQualities)
        fraction_5=sum(baseQualities[0:5])/totalbp*100
        fraction_10=sum(baseQualities[0:10])/totalbp*100
        fraction_15=sum(baseQualities[0:15])/totalbp*100
        fraction_20=sum(baseQualities[0:20])/totalbp*100
        boxplot_dic[sample] = [fraction_5,fraction_10,fraction_15,fraction_20]

    #draw boxplots
    sns.set_style("whitegrid")
    df=pd.DataFrame.from_dict(boxplot_dic)
    df=df.rename({0: "qual_5", 1: "qual_10", 2: "qual_15", 3: "qual_20"}, axis='index').transpose()
    sns.boxplot(data=df)
    plt.xlabel('Base qualities below given values')
    plt.ylabel('Percentage of base pairs (%)')
    plt.savefig("output/lowerQualityBPDistribution.png")
    plt.close()


def targetedAllMetricsProcess(infile):
    logging.info("Reading overall metrics file..")
    with open(infile,"r") as f:
        for line in f:
            if line.startswith("sample") and "bait_region" in line:
                header=line.rstrip().split("\t")


def main():

    parser = argparse.ArgumentParser(description='Script to plot useful data from a genome coverage analysis performed with GAT4K4. Matplotlib is required')
    parser.add_argument(dest='path',help='Path to the directory which holds all the output files.')
    parser.add_argument(dest='analysis_type', metavar='analysis_type', choices=("WGS", "targeted"), help="Type of the analysis performed. Choices: [WGS,targeted]")
    args = parser.parse_args()
    try:
        os.chdir(args.path)
        if not os.path.exists("output"):
            os.makedirs("output")
        if args.analysis_type == "WGS":
            wgs_all="final_colllectWgsMetrics_all.txt"
            if not os.path.exists(wgs_all):
                logging.error("Error: '{}' file doesn't exist in {} directory, which is required when WGS analysis is set.".format(wgs_all, args.path))
                sys.exit(1)
            elif len(glob.glob('*.histo')) == 0:
                logging.error("Error: No histogram files '*.histo' in {} directory. They are required in WGS experiments.".format(args.path))
                sys.exit(1)
            elif len(glob.glob('*.histo')) != (wccount(wgs_all) - 1):
                logging.error("Error. Number of '.histo' files must be the same as the number of samples described in the {} file.".format(wgs_all))
                sys.exit(1)
            #else:
            #    for f in glob.glob('*.histo'):

        elif args.analysis_type == "targeted":
            targeted_all = "final_collectHSmetrics_all.txt"
            perTarget_all = "final_perTargetCoverage.txt"
            if not os.path.exists(targeted_all) or not os.path.exists(perTarget_all):
                logging.error("Error: '{}' and '{}' must exist in {} directory. They are required when targeted analysis is set.".format(targeted_all, perTarget_all,args.path))
                sys.exit(1)
            elif countOcurrences(perTarget_all, "##") != (wccount(targeted_all) - 1):
                logging.error("Error. Number of samples differs in {} and {} files.".format(targeted_all, perTarget_all))
                sys.exit(1)
            else:
                processTargetedExperiment(targeted_all, perTarget_all)

    except OSError:
        logging.error("Error: {} is not a directory!".format(args.path))

if __name__ == "__main__":
    main()