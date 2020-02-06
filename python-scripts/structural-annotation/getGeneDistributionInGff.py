__author__ = 'pedro'

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
from collections import OrderedDict
from collections import Counter
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

def createGenomeDict(genome_table):
    base_dict=OrderedDict()

    with open(genome_table, 'r') as infile:
        for line in infile:
            base_dict[line.split()[0]] = line.rstrip().split()[1]

    print("Number of scaffolds in genome:\t%i" % len(base_dict))
    return base_dict

def processGff(gff_file):
    c = Counter()
    gene_length_strand = []
    with open(gff_file,'r') as infile:
        for line in infile:
            if "#" not in line and line.split()[2] == "gene":
                features= line.rstrip().split()
                c.update({features[0]: 1})
                gene_length_strand.append((features[8],int(features[4])-int(features[3]), features[6]))


    counter_strand = Counter(elem[2] for elem in gene_length_strand)
    avg_gene_length = round(sum(int(elem[1]) for elem in gene_length_strand)/sum(c.values()),2)
    longest_gene = max(int(elem[1]) for elem in gene_length_strand)
    shortest_gene = min(int(elem[1]) for elem in gene_length_strand)

    print("Number of scaffolds with genes predicted:\t%i" % len(list(c)))
    print("Number of genes predicted:\t%i" % sum(c.values()))
    print("Number of genes predicted in forward strand:\t%i" % counter_strand['+'])
    print("Number of genes predicted in reverse strand:\t%i" % counter_strand['-'])
    print("Average gene length:\t%i" % avg_gene_length)
    print("Maximum gene length\t%i" % longest_gene)
    print("Minimum gene length\t%i" % shortest_gene)

    return c, gene_length_strand

def writeOutFiles(base_dict, counter, geneLengthList,outbasename):

    lenght_vs_number = OrderedDict()
    with open(outbasename + "_individualScaffInfo.txt", 'w') as outfile:
        outfile.write("#scaffold_id\tlength\tpredicted_genes\n")
        for k,v in base_dict.items():
            if k in counter:
                outfile.write(k + "\t" + v + "\t" + str(counter[k]) + "\n")
                lenght_vs_number[int(v)] = counter[k]
            else:
                outfile.write(k + "\t" + v + "\t" + "0" + "\n")
                lenght_vs_number[int(v)] = 0
    outfile.close()

    with open(outbasename + "_individualGeneInfo.txt", 'w') as outfile:
        outfile.write("#gene_id\tlength\tstrand\n")
        for val in geneLengthList:
            outfile.write(val[0] + "\t" + str(val[1]) + "\t" + val[2] + "\n")
    outfile.close()

    return lenght_vs_number

def drawScatterPlot(lengthVSnumber, outbasename):

    outScatter = outbasename + "_scatterGenesInScaffolds.png"
    print("Drawing plot..")
    fig,ax = plt.subplots()
    for length,gene_number in lengthVSnumber.items():
        if gene_number == 0:
            ax.scatter(length,gene_number,color='red',s=25,alpha=0.6)
        else:
            ax.scatter(length,gene_number,color='grey',s=25,alpha=0.6)


    #draw subplot
    print("Drawing subplot..")
    axins = zoomed_inset_axes(ax, 7.5, loc=2) # zoom-factor: 2.5, location: upper-left
    for length,gene_number in lengthVSnumber.items():
        if gene_number == 0:
            gene_false = axins.scatter(length,gene_number,color='red',s=25,alpha=0.6)
        else:
            gene_true = axins.scatter(length,gene_number,color='grey',s=25,alpha=0.6)

    x1, x2, y1, y2 = 1000,100000 , -1, 30
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)

    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")


    plt.yticks(visible=False)
    plt.xticks(visible=False)
    axins.set_xlim(1000,100000)
    ax.set_xlabel("Scaffolds lenght")
    ax.set_ylabel("Number of predicted genes")
    ax.set_xlim(1000,max(lengthVSnumber.keys()) + 10000)
    ax.set_ylim(-1,max(lengthVSnumber.values()) + 40)
    ax.legend((gene_true,gene_false),('Scaffolds with genes predicted','Scaffolds without genes predicted'),
           scatterpoints=1,loc=1,ncol=1,fontsize=10)
    plt.savefig(outScatter)
    plt.close()
    print("Done")

def main():

    parser = argparse.ArgumentParser(description='Script to produce information about gene prediction distribution over a gff file.')
    parser.add_argument(dest='gff_file', metavar='gff', help='Annotation file to be analyzed.')
    parser.add_argument(dest='genome_table',metavar='genome',help='Genome table displaying scaffold and their length (one per line).')
    parser.add_argument(dest='output_basename',metavar='basename',help='Basename to write the output files.')
    args = parser.parse_args()

    base_dict = createGenomeDict(args.genome_table)
    counter, geneLengthList = processGff(args.gff_file)
    lengthVSnumber = writeOutFiles(base_dict,counter,geneLengthList,args.output_basename)
    drawScatterPlot(lengthVSnumber,args.output_basename)
if __name__ == "__main__":
    main()
