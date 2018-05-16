import argparse
from pybedtools import BedTool
import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO,  format='%(asctime)s %(message)s')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns

def addFeatureLengthToScoreCol(feature,col,posF,posS,delimeter):
    feature[4] = len(feature)
    if delimeter in feature[col-1]:
        feature[3] = feature[col-1].split(delimeter)[posF] + delimeter + feature[col-1].split(delimeter)[posS]
    else:
        print(delimeter)
        print("Delimiter provided not present in input bed.")
        exit(1)
    return feature

def filterByLength(feature,length):
    if len(feature) < length:
        return feature

def plotFeaturesLength(bedobj,columncol,posF,posS,delimeter,outbasename):
    df = bedobj.to_dataframe(names=['chrom', 'start', 'stop', 'name', 'length', 'strand'])
    df[['name', 'subfeature']] = df.name.str.split(delimeter, expand=True)
    df = df.groupby(['name', 'subfeature'],as_index = False)['length'].agg('sum')


    fig, (ax1, ax2) = plt.subplots(figsize=(10,8),ncols=1,nrows=2)
    sns.factorplot(x="name", y="length", hue="subfeature", data=df,
                       size=6, kind="bar", palette="muted",ax=ax1).despine(left=True)
    sns.factorplot(x="name", y="length", hue="subfeature", data=df,
                       size=6, kind="bar", palette="muted",ax=ax2).despine(left=True)


    ax2.set_ylim(0, 130000)
    #ax2.xaxis.tick_top()
    ax2.spines['top'].set_visible(False)
    ax2.xaxis.tick_bottom()
    ax2.set_xticklabels(df['name'].unique(), rotation=90)
    ax2.legend().set_visible(False)
    ax2.xaxis.label.set_visible(False)
    ax2.yaxis.label.set_visible(False)

    ax1.yaxis.label.set_visible(False)
    ax1.xaxis.label.set_visible(False)
    ax1.set_ylim(200000, df['length'].max() + 5000)
    ax1.spines['bottom'].set_visible(False)
    ax1.set_xticks([])
    # ax1.tick_params(labeltop='off')


    fig.text(0.04, 0.5, 'Length(bp)', va='center', rotation='vertical')
    d = .005  # how big to make the diagonal lines in axes coordinates
    # arguments to pass to plot, just so we don't keep repeating them
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((-d, +d), (-d, +d), **kwargs)  # top-left diagonal
    ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

    # What's cool about this is that now if we vary the distance between
    # ax and ax2 via f.subplots_adjust(hspace=...) or plt.subplot_tool(),
    # the diagonal lines will move accordingly, and stay right at the tips
    # of the spines they are 'breaking'
    fig.subplots_adjust(hspace=0.05)
    fig.savefig(outbasename + "_plot.pdf", format='pdf')
def processBed(inbed,col,posF,posS,delimeter,outbasename):
    logging.info("Creating BedTool object.")
    myBedTool = BedTool(inbed)
    logging.info("Retrieving feature length.")
    bedobj=myBedTool.each(addFeatureLengthToScoreCol,col,posF,posS,delimeter).saveas(outbasename + ".bed")
    plotFeaturesLength(bedobj,col,posF,posS,delimeter,outbasename)


def main():
    parser = argparse.ArgumentParser(description='Script to sum up the size of all subfeatures on a set of genes. Feature and subfeaures must coexist in the same string.'
                                                 'E.g ENSG00000178209.10_intron')
    parser.add_argument(dest='bed', help='Path to the bed file')
    parser.add_argument("-c", "--column", required=True, type=int,help='Column number in the bed file for which discriminative gene strings are represented')
    parser.add_argument('-d','--delimiter', required=True,help='Delimiter character to split feature/subfeature names')
    parser.add_argument("-p", "--positionFeature", required=True, type=int,help='Position index (0-based) where gene names are located in the "-c" columns, when splitted by "-d" delimiter')
    parser.add_argument("-s", "--positionSubfeature", required=True, type=int,help='Position index (0-based) where subfeatures names (e.g exon,intron) are located in the "-c" column, when splitted by "-d" delimiter')
    parser.add_argument("-o", "--output", required=True, help='Basename to write output files.')

    args = parser.parse_args()

    processBed(args.bed, args.column,args.positionFeature,args.positionSubfeature,args.delimiter,args.output)

if __name__ == "__main__":
    main()
