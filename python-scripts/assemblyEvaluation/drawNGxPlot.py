from Bio import SeqIO
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab
import matplotlib.cm as cmx
import matplotlib.colors as colors
import sys
import os
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
from collections import OrderedDict
from collections import defaultdict

def declareGLobal():
    global dict_glb
    dict_glb=OrderedDict()

def processFasta(fasta):
    logging.info("Processing %s file.." % fasta)
    handle = open(fasta,'rU')
    sequences = SeqIO.parse(handle,'fasta')
    basename= os.path.basename(fasta)
    assemblyname=""
    if not '.fasta' or not '.fna' or not '.fa' in basename:
        logging.info("You sure you provided fasta files? Can't find any of the required extensions.")
        exit(1)
    elif '.fasta' in basename:
        assemblyname=basename.split(".fasta")[0]
    elif '.fna' in basename:
        assemblyname=basename.split(".fna")[0]
    elif '.fa' in basename:
        assemblyname=basename.split(".fa")[0]

    gntable={}
    logging.info("\tCalculating lengths..")
    for record in sequences:
        if record.description in gntable:
            logging.error("Repetitive sequence header (%s) in %s file ..\nExiting." % (record.description,basename))
            exit(1)
        else:
            gntable[record.description] = len(record.seq)

    handle.close()
    sequences.close()
    if basename in dict_glb:
        logging.error("You provided the %s assembly more than once. Please remove one of them from the args list" % basename)
        exit(1)
    else:
        dict_glb[assemblyname] = gntable

def getCumLength(genomeSize,stepSize):
    logging.info("Calculating metrics.")
    genomeFractions=np.arange(0,101,stepSize)
    ngx_dict=defaultdict(list)
    for assembly,contigLen in dict_glb.items():
        sort_list=sorted(contigLen.values(),reverse=True)
        for percentage in genomeFractions:
            ngx=getNGx(genomeSize,percentage,sort_list)
            ngx_dict[assembly].append(ngx)

    return ngx_dict,genomeFractions


def drawplot(finaldict,genomeFractions):
    logging.info("Drawing plot..")
    fig, ax = plt.subplots()
    lists=sorted(finaldict.items())
    assemblienames,y = zip(*lists)

    #x labels
    fraction2NG=["NG"+str(p) for p in genomeFractions]
    with open('NGtable.txt','w') as outfile:
        outfile.write("#Name\t" + '\t'.join(fraction2NG)+"\n")
        for assembly, ngx in finaldict.items():
            outfile.write(assembly +"\t" + "\t".join(str(j) for j in ngx)+ "\n")
    bins=7
    if len(fraction2NG) >= bins:
        arr=np.array_split(fraction2NG,bins)
        ticks=[]
        for a in arr:
            ticks.append(a[0])
    else:
        bins=len(fraction2NG)
        ticks=fraction2NG


    #y labels
    #colors
    cmap=get_cmap(len(assemblienames))


    #draw plot
    for i in range(0,len(y)):
        col=cmap(i)
        plt.xticks(genomeFractions,ticks)

        lines=plt.plot(genomeFractions, [j for j in y[i]], label=assemblienames[i])
        plt.setp(lines, color=col, linewidth=1.0)

    handles, labels = ax.get_legend_handles_labels()
    lgd=ax.legend(handles, labels,loc='upper center', ncol=2,bbox_to_anchor=(0.5,-0.1) )

    ax.locator_params(nbins=7, axis='x')
    ax.set_xlabel('NG fraction')
    ax.set_ylabel('Contig size (bp)')
    plt.axvline(x=50,linestyle='dashed',color='darkgrey')
    for legobj in lgd.legendHandles:
        legobj.set_linewidth(2.0)
    fig.savefig('NGxplot.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close()
    logging.info("Done")


def get_cmap(N):
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='PiYG')
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color


def getNGx(genomeSize, percentage,sortlist):
    s = genomeSize
    limit = int(genomeSize * (100.0 - percentage) / 100)
    for len in sortlist:
        s -= len
        if s <= int(limit):
            ngx = len
            return ngx

    return None

def main():
    parser = argparse.ArgumentParser(description='Script to draw cumulative assembly plots from fasta(s) file(s).')
    parser.add_argument(dest='assemblies', nargs='+', metavar='assembly_fasta', help='List of fasta files to draw the plot. Script will look for any fasta extension [fasta,fa,fna] to split and '
                                                                                     'assign each assembly name.')
    parser.add_argument('-g',metavar='--genomeSize', required=True, type=int, help='Genome size value to perform the analysis')
    parser.add_argument('-l', '--list', action='store_true', help='Input is a file listing all the fasta together, one per line.')
    parser.add_argument('-f', metavar='--fractions', default=5, type=int,help='Step size to calculate NG[x]. Default: 5.')
    args = parser.parse_args()

    declareGLobal()
    if len(args.assemblies) > 1 and args.list:
        logging.info("When '-l' set, please don't provide more than 1 file in the positional arguments." )
        exit(1)

    elif args.list:
        with open(args.assemblies[0], 'r') as listFasta:
            for line in listFasta:
                l=line.rstrip()
                processFasta(l)
            final,genomeFractions=getCumLength(args.g,args.f)
            drawplot(final,genomeFractions)

    else:
        for file in args.assemblies:
            processFasta(file)
        final,genomeFractions=getCumLength(args.g,args.f)
        drawplot(final,genomeFractions)

if __name__ == "__main__":
    main()