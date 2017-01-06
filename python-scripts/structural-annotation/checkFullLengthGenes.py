import argparse
import logging
import sys
import collections
import operator
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')

def processGff(gff):

    dict_utr = {}
    dict_gene_coord_strand = {}
    dict_numberOfStartStop = {}
    dict_gene_transcript ={}
    nostart,nostop,notstartstop,oneone, several=[],[],[],[],[]
    firstStart=False
    firstStop=False
    existStart=True
    existStop=True
    i=0
    with open(gff, 'r') as infile:
        previou_g_id = ""
        for line in infile:
            if not line.startswith("#"):
                if line.split("\t")[2] == "gene":

                    if previou_g_id in dict_numberOfStartStop.keys():
                        if existStart == False and existStop == False:
                            dict_utr[previou_g_id][0] = "no_startCodon"
                            dict_utr[previou_g_id][1] = "no_stopCodon"
                        elif existStart == False:
                            dict_utr[previou_g_id][0] = "no_startCodon"
                        elif existStop == False:
                            dict_utr[previou_g_id][1] = "no_stopCodon"

                        if dict_numberOfStartStop[previou_g_id][0] == 0 and dict_numberOfStartStop[previou_g_id][1] == 0:
                            logging.info("IMPORTANT: %g gene does not have both start and stop codons")
                            notstartstop.append(previou_g_id)

                        elif dict_numberOfStartStop[previou_g_id][0] == 0:
                            nostart.append(previou_g_id)

                        elif dict_numberOfStartStop[previou_g_id][1] == 0:
                            nostop.append(previou_g_id)

                        elif dict_numberOfStartStop[previou_g_id][0] == 1 and dict_numberOfStartStop[previou_g_id][1] == 1:
                            oneone.append(previou_g_id)

                        elif dict_numberOfStartStop[previou_g_id][0] > 1 or dict_numberOfStartStop[previou_g_id][1] > 1:
                            several.append(previou_g_id)



                    features = line.rstrip().split("\t")
                    dict_gene_coord_strand[features[8]] = [features[3],features[4],features[6]]
                    dict_utr[features[8]] = ["",""]
                    dict_numberOfStartStop[features[8]] = [0,0]
                    dict_gene_transcript[features[8]] = []
                    previou_g_id = features[8]
                    firstStart=True
                    firstStop=True
                    existStart=False
                    existStop=False


                elif "codon" in line.split("\t")[2]:
                    feature = line.split("\t")
                    attr_list = line.split("\t")[8]
                    gene = [g for g in attr_list.split(";") if "gene_id" in g]
                    if len(gene) == 0:
                        logging.error("No gene_id attribute in the gff")
                        exit(1)
                    g_id = gene[0].split("\"")[1]

                    if feature[2] == "start_codon":
                        existStart=True
                        dict_numberOfStartStop[g_id][0] += 1
                        #forward
                        if firstStart == True and dict_gene_coord_strand[g_id][2] == "+":

                            if int(feature[3]) == int(dict_gene_coord_strand[g_id][0]):
                                dict_utr[g_id][0] = "no_5primeUTR"

                            elif int(feature[3]) > int(dict_gene_coord_strand[g_id][0]):
                                length= int(feature[3]) - int(dict_gene_coord_strand[g_id][0])
                                dict_utr[g_id][0] = "5primerUTR_" + str(length)

                            firstStart=False

                        elif firstStart == False and dict_gene_coord_strand[g_id][2] == "+":

                            logging.info("%s has multiple start codon.." % g_id)
                            if int(feature[3]) == int(dict_gene_coord_strand[g_id][0]):
                                dict_utr[g_id][0] = "no_5primeUTR"

                        #reverse
                        elif firstStart==True and dict_gene_coord_strand[g_id][2] == "-":

                            if feature[4] == dict_gene_coord_strand[g_id][1]:
                                dict_utr[g_id][0] = "no_5primeUTR"
                                firstStart=False

                            elif int(feature[4]) < int(dict_gene_coord_strand[g_id][1]):
                                length= int(dict_gene_coord_strand[g_id][1]) - int(feature[4])
                                dict_utr[g_id][0] = "5primerUTR_" + str(length)
                                firstStart=False

                        elif firstStart==False and dict_gene_coord_strand[g_id][2] == "-":#when genes in negative strand we need to loo until last start codon appeared
                            logging.info("%s has multiple start codon.." % g_id)
                            if int(feature[4]) < int(dict_gene_coord_strand[g_id][1]):
                                length= int(dict_gene_coord_strand[g_id][1]) - int(feature[4])
                                dict_utr[g_id][0] = "5primerUTR_" + str(length)

                            elif int(feature[4]) == int(dict_gene_coord_strand[g_id][1]):
                                dict_utr[g_id][0] = "no_5primeUTR"

                    ###############################################
                    elif feature[2] == "stop_codon":
                        existStop=True
                        dict_numberOfStartStop[g_id][1] += 1
                        #forward
                        if firstStop==True and dict_gene_coord_strand[g_id][2] == "+":

                            if int(feature[4]) == int(dict_gene_coord_strand[g_id][1]):
                                dict_utr[g_id][1] = "no_3primeUTR"

                            elif int(feature[4]) < int(dict_gene_coord_strand[g_id][1]):
                                length= int(dict_gene_coord_strand[g_id][0]) - int(feature[4])
                                dict_utr[g_id][1] = "3primerUTR_" + str(length)
                            firstStop=False

                        elif firstStop==False and dict_gene_coord_strand[g_id][2] == "+":#when genes in positive strand we need to loop until last stop codon appears
                            logging.info("%s has multiple stop codon.." % g_id)
                            if int(feature[4]) < int(dict_gene_coord_strand[g_id][1]):
                                length= int(dict_gene_coord_strand[g_id][1]) - int(feature[4])
                                dict_utr[g_id][1] = "3primerUTR_" + str(length)
                            elif int(feature[4]) == int(dict_gene_coord_strand[g_id][1]):
                                dict_utr[g_id][1] = "no_3primeUTR"
                        #reverse
                        elif firstStop==True and dict_gene_coord_strand[g_id][2] == "-":

                            if feature[3] == dict_gene_coord_strand[g_id][0]:
                                dict_utr[g_id][1] = "no_3primeUTR"
                                firstStop=False

                            elif int(feature[3]) > int(dict_gene_coord_strand[g_id][1]):
                                length=  int(feature[3]) - int(dict_gene_coord_strand[g_id][1])
                                dict_utr[g_id][1] = "3primerUTR_" + str(length)
                                firstStart=False

                        elif firstStop == False and dict_gene_coord_strand[g_id][2] == "-":
                            logging.info("%s has multiple stop codon.." % g_id)
                            if feature[3] == dict_gene_coord_strand[g_id][0]:
                                dict_utr[g_id][1] = "no_3primeUTR"

                elif line.split("\t")[2] == "transcript":

                    dict_gene_transcript[previou_g_id].append(line.split("\t")[8].rstrip())


        #last gene
        if previou_g_id in dict_numberOfStartStop.keys():
            if existStart == False and existStop == False:
                dict_utr[previou_g_id][0] = "no_startCodon"
                dict_utr[previou_g_id][1] = "no_stopCodon"
            elif existStart == False:
                dict_utr[previou_g_id][0] = "no_startCodon"
            elif existStop == False:
                dict_utr[previou_g_id][1] = "no_stopCodon"

            if dict_numberOfStartStop[previou_g_id][0] == 0 and dict_numberOfStartStop[previou_g_id][1] == 0:
                logging.info("IMPORTANT: %s gene does not have both start and stop codons" % previou_g_id)
                notstartstop.append(previou_g_id)

            elif dict_numberOfStartStop[previou_g_id][0] == 0:
                nostart.append(previou_g_id)

            elif dict_numberOfStartStop[previou_g_id][1] == 0.:
                nostop.append(previou_g_id)

            elif dict_numberOfStartStop[previou_g_id][0] == 1 and dict_numberOfStartStop[previou_g_id][1] == 1:
                oneone.append(previou_g_id)

            elif dict_numberOfStartStop[previou_g_id][0] > 1 or dict_numberOfStartStop[previou_g_id][1] > 1:
                several.append(previou_g_id)


        logging.info("\n################################\nNumber of genes without start codon\t%i" %len(nostart))
        logging.info("Number of genes without stop codon\t%i" %len(nostop))
        logging.info("Number of genes without start and stop codon\t%i\n" %len(notstartstop))
        logging.info("Number of genes with exactly one start and one stop codon\t%i" %len(oneone))
        logging.info("Number of genes with several start or stop codons\t%i\n" %len(several))


        noutr5,noutr3,utr5prime,utr3prime,noStart,noStop=0,0,0,0,0,0

        logging.info("####UTR analysis####")
        for k,v in dict_utr.items():
            if "no_5prime" in v[0]:
                noutr5+=1
            elif "no_start" in v[0]:
                noStart+=1
            else:
                utr5prime+=1
                #print(k)
            if "no_3prime" in v[1]:
                noutr3+=1

            elif "no_stop" in v[1]:
                noStop+=1
            else:
                utr3prime+=1
                #print(k)

        logging.info("Number of genes with 5 prime UTR\t%i" % utr5prime)
        logging.info("Number of genes with 3 prime UTR\t%i" % utr3prime)
        logging.info("Genes with no 5 prime\t%i" % noutr5)
        logging.info("Genes with no 3 prime\t%i" % noutr3)
        logging.info("Genes with unknown 5prime UTR prediction (no start codon)\t%i" % noStart)
        logging.info("Genes with unknown 3prime UTR prediction (no stop codon)\t%i" % noStop)

        with open("gene_information.txt", 'w') as outfile:
            outfile.write("#gene_id\tstart_codon\tstop_codon\tnumber_start\tnumber_stop\tnumber_transcripts\ttranscripts_name\n")
            for k,v in sorted(dict_gene_transcript.items(), key = lambda item : len(item[1]), reverse=True):
                start=["no" if k in nostart or k in notstartstop else "yes"]
                stop=["no" if k in nostop or k in notstartstop else "yes"]
                outfile.write(k + "\t" + start[0] + "\t" + stop[0] + "\t" + str(dict_numberOfStartStop[k][0]) + "\t" + str(dict_numberOfStartStop[k][1])
                              + "\t" + str(len(v)) + "\t" + ';'.join(v) + "\n")

        for k,v in dict_gene_transcript.items():
            i+=len(dict_gene_transcript[k])
        #print(i)
        Count = collections.Counter([len(v) for k,v in dict_gene_transcript.items()])
        print(Count)
def main():
    parser = argparse.ArgumentParser(description='Script to look for start, stop and UTR region in gff3 files.')
    parser.add_argument(dest='gff_file', metavar='gff', help='Annotation file to be analyzed.')
    #parser.add_argument('-w', '--write', help='Write new gff file with UTR regions added.')
    args = parser.parse_args()

    processGff(args.gff_file)

if __name__ == "__main__":
    main()

