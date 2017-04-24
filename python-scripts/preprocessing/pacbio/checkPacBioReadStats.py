import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import SeqIO
from collections import OrderedDict
import os
import numpy as np
import seaborn as sns
import pandas as pd
import logging
import sys
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')

def declareGLobal():
    global dataStruct
    dataStruct=OrderedDict()

def processFastaFiles(inputFile):

    smrtcell_id = os.path.basename(inputFile).split("_")[3]

    if smrtcell_id in dataStruct:
        logging.info("Processing another file of existing %s smrtcell" % smrtcell_id)
        indivSmrtcellDic = dataStruct[smrtcell_id]
    else:
        logging.info("Processing first file of %s smrtcell" % smrtcell_id)
        indivSmrtcellDic = OrderedDict()
    rq = 0

    handle = open(inputFile,'rU')
    sequences = SeqIO.parse(handle,'fasta')
    for record in sequences:

        hole_number=record.description.split("/")[1]
        if not "RQ=" in record.description:
            logging.info("read quality tag may not be present in fasta headers. [%s]" % record.description)
            exit(1)
        else:
            rq=record.description.split("RQ=")[1]

        if not hole_number in indivSmrtcellDic:
            indivSmrtcellDic[hole_number] = [rq]

        indivSmrtcellDic[hole_number].append(len(record.seq))


    dataStruct[smrtcell_id] = indivSmrtcellDic
    handle.close()

def generateStats(outbasename):
    logging.info("Generating stats..")
    zmw_perSmrtcell,subread_quality=[],[]
    throuhtput_persmrtcell,singlepassSubreads,passes_perZMW= {},{},{}


    with open(outbasename + "_stats.txt", 'w') as out1:
        out1.write("Total number of smrt cells:\t%i\n" % len(dataStruct))
        for smrtcell, dict_zmw in dataStruct.items():
            zmw_perSmrtcell.append(len(dict_zmw))
            throughput=0

            for zmw in dict_zmw:
                id = smrtcell + "_" + zmw
                passes_perZMW[id] = (len(dict_zmw[zmw]) - 1)
                if len(dict_zmw[zmw]) - 1 == 1:
                    singlepassSubreads[id] = dict_zmw[zmw][1]

                for subread in dict_zmw[zmw][1:]:

                    subread_quality.append((int(subread),float(dict_zmw[zmw][0])))
                    throughput+=subread

            throuhtput_persmrtcell[smrtcell] = throughput

        out1.write("Average number of ZMW per smrtcell:\t%.2f\n" % np.mean(zmw_perSmrtcell))
        out1.write("Sequence throughput (Gb):\t%.4f\n" % (sum(throuhtput_persmrtcell.values())/1000000000))
        out1.write("Total number of subreads:\t%i\n" % len(subread_quality))
        out1.write("Average subread length (bp):\t%.2f\n" % (sum([ln[0] for ln in subread_quality])/len(subread_quality)))
        out1.write("Average subread quality:\t%.4f\n" % (sum(rq[1] for rq in subread_quality)/len(subread_quality)))
        out1.write("Average number of passes represented in subreads:\t%.2f\n" % np.mean(list(passes_perZMW.values())))
        out1.write("Number of subreads coming from single pass ZMW:\t%i\n" % len(singlepassSubreads))
        out1.write("Maximum number of passes observed in a polymerase read:\t%i [%s]\n" % (np.max(list(passes_perZMW.values())),max(passes_perZMW, key=passes_perZMW.get)))


        len_list= np.asarray([int(i[0]) for i in subread_quality])
        hist,bin_edges = np.histogram(len_list,50,(0,100000))

        out1.write("Number of reads per bin:\n")
        for i in range(0,len(hist)):
            out1.write("n >= %i\t%i\n" % (bin_edges[i],np.sum(list((hist[i:])))))

    out1.close()
    with(open(outbasename + "_smrtcellsThroughput.tsv",'w')) as out2:
        out2.write("#smrtcell_id\tthroughput(Gb)\n")
        for k,v in throuhtput_persmrtcell.items():
            out2.write("%s\t%.4f\n" % (k,v/1000000000))

    logging.info("Drawing plots..")
    ###subread length hist
    logging.info("\tHistogram..")
    lengths=[i[0] for i in subread_quality]
    plt.hist(lengths, bins= 200, color='#cdb79e', range=[0,50000],histtype='stepfilled')
    plt.xlabel('Read lenght (bp)')
    plt.ylabel('Counts')
    #plt.axvline(np.mean(lengths), color='#838b83', linestyle='dashed', linewidth=1)
    #plt.axvline(np.median(lengths), color='#cd853f', linestyle='dashed', linewidth=1)
    plt.savefig(outbasename + "_readLenDistribution.png")
    plt.close()


    logging.info("\tScatter plot..")
    #scatter plot
    df = pd.DataFrame(subread_quality,columns=("subread_length","subread_quality"))
    if df.shape[0] > 500000:
        logging.info("Number of subreads is higher than 500,000, random subsampling will be performed.")
        df = df.sample(n=500000)
    sns.jointplot(x='subread_length', y='subread_quality',data=df,kind='kde')
    plt.savefig(outbasename + "_scatterQualities.png")
    plt.close()

    logging.info("\tSingle pass density plot..")
    #single pass subreads plot
    single_subread_len=list(singlepassSubreads.values())
    sns.kdeplot(np.array(single_subread_len), bw=0.5)
    plt.title("Density of single pass subreads length")
    plt.savefig(outbasename + "_single.png")
    plt.close()



def main():

    parser = argparse.ArgumentParser(description='Script to check overall stats from pacbio subreads data in fasta format.')
    parser.add_argument(dest='fasta_file', metavar='fasta', nargs="+", help='Fasta files to be processed.')
    parser.add_argument('-o', metavar = 'output_basename',required=True,help='Basename to write the output files.')
    parser.add_argument('-l', '--list', action='store_true', help='Input is a file listing all the fasta together, one per line.')
    args = parser.parse_args()

    declareGLobal()
    if len(args.fasta_file) > 1 and args.list:
        logging.info("When '-l' set, please don't provide more than 1 file in the positional arguments." )
        exit(1)

    elif args.list:
        with open(args.fasta_file[0], 'r') as listFasta:
            for line in listFasta:
                l=line.rstrip()
                processFastaFiles(l)
            generateStats(args.o)
    else:
        for file in args.fasta_file:
            processFastaFiles(file)
        generateStats(args.o)
if __name__ == "__main__":
    main()