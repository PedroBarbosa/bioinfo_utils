import argparse
import logging
import sys
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
    endstartUTRpred=False
    endstopUTRpred=False
    with open(gff, 'r') as infile:
        previou_g_id = ""
        for line in infile:
            if not line.startswith("#"):
                if line.split("\t")[2] == "transcript":

                    if previou_g_id in dict_numberOfStartStop.keys():
                        if existStart == False and existStop == False:
                            dict_utr[previou_g_id][0] = "no_startCodon"
                            dict_utr[previou_g_id][1] = "no_stopCodon"
                        elif existStart == False:
                            dict_utr[previou_g_id][0] = "no_startCodon"
                        elif existStop == False:
                            dict_utr[previou_g_id][1] = "no_stopCodon"

                        if dict_numberOfStartStop[previou_g_id][0] == 0 and dict_numberOfStartStop[previou_g_id][1] == 0:
                            logging.info("IMPORTANT: %s transcript does not have both start and stop codons" % previou_g_id)
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
                    endstartUTRpred=False
                    endstopUTRpred=False

                elif "codon" in line.split("\t")[2]:
                    feature = line.split("\t")
                    attr_list = line.split("\t")[8]
                    gene = [g for g in attr_list.split(";") if "transcript_id" in g]
                    if len(gene) == 0:
                        logging.error("No transcript_id attribute in the gff")
                        exit(1)
                    g_id = gene[0].split("\"")[1]

                    if feature[2] == "start_codon":
                        existStart=True
                        dict_numberOfStartStop[g_id][0] += 1
                        #forward
                        if firstStart == True and dict_gene_coord_strand[g_id][2] == "+" and endstartUTRpred==False:

                            if int(feature[3]) == int(dict_gene_coord_strand[g_id][0]):
                                dict_utr[g_id][0] = "no_5primeUTR"
                                endstartUTRpred=True
                            elif int(feature[3]) > int(dict_gene_coord_strand[g_id][0]):
                                length= int(feature[3]) - int(dict_gene_coord_strand[g_id][0])
                                dict_utr[g_id][0] = "5primerUTR_" + str(length)

                            firstStart=False

                        elif firstStart == False and dict_gene_coord_strand[g_id][2] == "+" and endstartUTRpred==False:

                            logging.info("%s has multiple start codon.." % g_id)
                            if int(feature[3]) == int(dict_gene_coord_strand[g_id][0]):
                                dict_utr[g_id][0] = "no_5primeUTR"
                                endstartUTRpred=True
                        #reverse
                        elif firstStart==True and dict_gene_coord_strand[g_id][2] == "-" and endstartUTRpred==False:

                            if feature[4] == dict_gene_coord_strand[g_id][1]:
                                dict_utr[g_id][0] = "no_5primeUTR"
                                firstStart=False
                                endstartUTRpred=True
                            elif int(feature[4]) < int(dict_gene_coord_strand[g_id][1]):
                                length= int(dict_gene_coord_strand[g_id][1]) - int(feature[4])
                                dict_utr[g_id][0] = "5primerUTR_" + str(length)
                                firstStart=False

                        elif firstStart==False and dict_gene_coord_strand[g_id][2] == "-" and endstartUTRpred==False:#when genes in negative strand we need to loo until last start codon appeared
                            logging.info("%s has multiple start codon.." % g_id)
                            if int(feature[4]) < int(dict_gene_coord_strand[g_id][1]):
                                length= int(dict_gene_coord_strand[g_id][1]) - int(feature[4])
                                dict_utr[g_id][0] = "5primerUTR_" + str(length)

                            elif int(feature[4]) == int(dict_gene_coord_strand[g_id][1]):
                                dict_utr[g_id][0] = "no_5primeUTR"
                                endstartUTRpred=True

                    ###############################################
                    elif feature[2] == "stop_codon":
                        existStop=True
                        dict_numberOfStartStop[g_id][1] += 1
                        #forward
                        if firstStop==True and dict_gene_coord_strand[g_id][2] == "+" and endstopUTRpred == False:

                            if int(feature[4]) == int(dict_gene_coord_strand[g_id][1]):
                                dict_utr[g_id][1] = "no_3primeUTR"
                                endstopUTRpred = True
                            elif int(feature[4]) < int(dict_gene_coord_strand[g_id][1]):
                                length= int(dict_gene_coord_strand[g_id][1]) - int(feature[4])
                                dict_utr[g_id][1] = "3primerUTR_" + str(length)
                            firstStop=False

                        elif firstStop==False and dict_gene_coord_strand[g_id][2] == "+" and endstopUTRpred == False:#when genes in positive strand we need to loop until last stop codon appears
                            logging.info("%s has multiple stop codon.." % g_id)
                            if int(feature[4]) < int(dict_gene_coord_strand[g_id][1]):
                                length= int(dict_gene_coord_strand[g_id][1]) - int(feature[4])
                                dict_utr[g_id][1] = "3primerUTR_" + str(length)
                            elif int(feature[4]) == int(dict_gene_coord_strand[g_id][1]):
                                dict_utr[g_id][1] = "no_3primeUTR"
                                endstopUTRpred = True
                        #reverse
                        elif firstStop==True and dict_gene_coord_strand[g_id][2] == "-" and endstopUTRpred == False:

                            if feature[3] == dict_gene_coord_strand[g_id][0]:
                                dict_utr[g_id][1] = "no_3primeUTR"
                                firstStop=False
                                endstopUTRpred=True

                            elif int(feature[3]) > int(dict_gene_coord_strand[g_id][0]):
                                length=  int(feature[3]) - int(dict_gene_coord_strand[g_id][0])
                                dict_utr[g_id][1] = "3primerUTR_" + str(length)
                                firstStart=False

                        elif firstStop == False and dict_gene_coord_strand[g_id][2] == "-" and endstopUTRpred == False:
                            logging.info("%s has multiple stop codon.." % g_id)
                            if feature[3] == dict_gene_coord_strand[g_id][0]:
                                dict_utr[g_id][1] = "no_3primeUTR"
                                endstopUTRpred = True

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
                logging.info("IMPORTANT: %s transcript does not have both start and stop codons" % previou_g_id)
                notstartstop.append(previou_g_id)

            elif dict_numberOfStartStop[previou_g_id][0] == 0:
                nostart.append(previou_g_id)

            elif dict_numberOfStartStop[previou_g_id][1] == 0.:
                nostop.append(previou_g_id)

            elif dict_numberOfStartStop[previou_g_id][0] == 1 and dict_numberOfStartStop[previou_g_id][1] == 1:
                oneone.append(previou_g_id)

            elif dict_numberOfStartStop[previou_g_id][0] > 1 or dict_numberOfStartStop[previou_g_id][1] > 1:
                several.append(previou_g_id)


        logging.info("\n################################\n")
        logging.info("Number of transcripts without start codon\t%i" %len(nostart))
        logging.info("Number of transcripts without stop codon\t%i" %len(nostop))
        logging.info("Number of transcripts without start and stop codon\t%i\n" %len(notstartstop))
        logging.info("Number of transcripts with exactly one start and one stop codon\t%i" %len(oneone))
        logging.info("Number of transcripts with several start or stop codons\t%i\n" %len(several))


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

        logging.info("Number of transcripts with 5 prime UTR (the existing start codons from these genes are distinct from the begin coordinate of it. Start codons of some transcripts of gene may be missing).\t%i" % utr5prime)
        logging.info("Number of transcripts with 3 prime UTR (the existing stop codons from these genes are distinct from the end coordinate of it. Stop codons of some transcripts of gene may be missing).\t%i" % utr3prime)
        logging.info("Transcripts with no 5 prime\t%i" % noutr5)
        logging.info("Transcripts with no 3 prime\t%i" % noutr3)
        logging.info("Transcripts with unknown 5prime UTR prediction (no start codon)\t%i" % noStart)
        logging.info("Transcripts with unknown 3prime UTR prediction (no stop codon)\t%i\n" % noStop)

        logging.info("Writing general output file..")
        with open("transcripts_information.txt", 'w') as outfile:
            outfile.write("#gene_id\tstart_codon\tstop_codon\tnumber_start\tnumber_stop\tnumber_transcripts\ttranscripts_name\t5primeUTR\t3primeUTR\n")
            for k,v in sorted(dict_gene_transcript.items(), key = lambda item : len(item[1]), reverse=True):
                start=["no" if k in nostart or k in notstartstop else "yes"]
                stop=["no" if k in nostop or k in notstartstop else "yes"]
                utr5=["no" if "no_" in dict_utr[k][0] else "yes_" + dict_utr[k][0].split("_")[1]]
                utr3=["no" if "no_" in dict_utr[k][1] else "yes_" + dict_utr[k][1].split("_")[1]]
                outfile.write(k + "\t" + start[0] + "\t" + stop[0] + "\t" + str(dict_numberOfStartStop[k][0]) + "\t" + str(dict_numberOfStartStop[k][1])
                              + "\t" + str(len(v)) + "\t" + ';'.join(v) + "\t" + utr5[0] + "\t" + utr3[0] + "\n")

        outfile.close()

        logging.info("Writing potential transcripts to remove from annotation (no stop, no start, not stop and start)..")
        with open("incomplete_transcripts.txt", "w") as outfile2:
            if notstartstop:
                outfile2.write("\n".join(notstartstop) + "\n")
            if nostart:
                outfile2.write("\n".join(nostart) + "\n")
            if nostop:
                outfile2.write("\n".join(nostop) + "\n")

def main():
    parser = argparse.ArgumentParser(description='Script to look for start, stop and UTR region in gff3 files.')
    parser.add_argument(dest='gff_file', metavar='gff', help='Annotation file to be analyzed.')
    #parser.add_argument('-w', '--write', help='Write new gff file with UTR regions added.')
    args = parser.parse_args()

    processGff(args.gff_file)

if __name__ == "__main__":
    main()

