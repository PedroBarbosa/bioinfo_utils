import argparse
import logging
import sys
import os
import numpy as np
from collections import defaultdict
from collections import OrderedDict
from collections import Counter
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')


class validateAnnotation():
    def __init__(self, software, outputprefix, genome_table, splitstring):
        self.software = software
        self.outputprefix = outputprefix
        self.genome_table = genome_table
        self.splitstring = splitstring

        self.scaffolds = set() #set of scaffolds with genes in the reference annotation
        self.genes = set() #set of genes in the reference annotation
        self.transcripts = set() # set of transcripts in the reference annotation
        self.allpresentscaffolds = set() #set of scaffolds with genes predicted in the set of assembled transcripts
        self.abundance = defaultdict(list)  # dictionary of all the transcripts TPM abundance estimation
        self.geneInterception = defaultdict(list) #dictionary of all present/absent genes
        self.individualGeneMatchingRate = OrderedDict() #list representing average of accurate genes predicted
        self.individualTranscriptMatchingRate = OrderedDict()  #list representing average of accurate transcripts predicted
        self.genesFullyCovered = set() #set of genes fully covered by RNA-seq reads (availalbe only when cov ref files are provided)
        self.transcriptsFullyCovered = set() #set of trasncripts fully covered by RNA-seq reads (availalbe only when cov ref files are provided)
        self.setStrands = set()
        self.indexDict = OrderedDict()
    def processReferenceAnnotGff(self,gff):
        logging.info("Processing reference annotation file..")
        with open(gff,'r') as infile:
            for line in infile:
                tmpline = line.rstrip().split("\t")
                if "#" not in line and tmpline[2] == "gene":
                    self.genes.add(tmpline[-1])
                    self.scaffolds.add(tmpline[0])

                elif "#" not in line and tmpline[2] == "transcript":
                    attr = tmpline[-1].split(";")
                    self.transcripts.add(attr[0])

    def processFromStringtie(self,gtf, index, numbfiles):
        logging.info("Reading stringtie transcripts file (%s).." % gtf)
        indivscafollds = set()
        indivgenes = set()
        indivnewgenes = set()
        indivtranscript = set()
        indivnewtranscript = set()
        with open(gtf, 'r') as infile:
            self.indexDict[index] = os.path.basename(gtf).split(self.splitstring)[0]
            logging.info("Updating global variables..")
            for line in infile:
                if not line.startswith("#"):
                    if line.split("\t")[2] == "transcript" or line.split("\t")[2] == "mRNA":

                        gene,new_gene,transcript,new_transcript,g_id,t_id="","","","","",""
                        tmpline = line.rstrip().split("\t")
                        scaffold_id = tmpline[0]
                        strand = tmpline[6]
                        attr = tmpline[8]
                        gene = [g for g in attr.split(";") if "ref_gene_id" in g]
                        transcript = [t for t in attr.split(";") if "reference_id" in t]
                        if not gene and not transcript:
                            new_gene = [g for g in attr.split(";") if "gene_id" in g]
                            new_transcript = [g for g in attr.split(";") if "transcript_id" in g]
                            g_id = new_gene[0].split("\"")[1]
                            t_id = new_transcript[0].split("\"")[1]
                        elif not gene and transcript or gene and not transcript:
                            logging.error("Error in Stringite output. Line %s should have either both reference_id and ref_gene_id or none." % line.rstrip())
                            exit(1)
                        else:
                            g_id = gene[0].split("\"")[1]
                            t_id = transcript[0].split("\"")[1]


                        self.allpresentscaffolds.add(scaffold_id)
                        indivscafollds.add(scaffold_id)
                        if gene:
                            indivgenes.add(g_id)
                            indivtranscript.add(t_id)
                            self.setStrands.add((g_id,strand))
                            tpm = [g for g in attr.split(";") if "TPM" in g]
                            if not self.abundance[t_id]:
                                self.abundance[t_id] = [0] * numbfiles
                            if not self.geneInterception[g_id]:
                                self.geneInterception[g_id] = ["no"] * numbfiles
                            self.abundance[t_id][index] = tpm[0].split("\"")[1]
                            self.geneInterception[g_id][index] = "yes"
                        else:
                            indivnewgenes.add(g_id)
                            indivnewtranscript.add(t_id)

            logging.info("Done.")
            self.writeIndividualStats(gtf,indivscafollds,indivgenes,indivnewgenes,indivtranscript,indivnewtranscript)


    def processCovRefs(self,covfile):
        fully_cov_genes = set()
        fully_cov_transcripts = set()
        logging.info("Processing covRef file..")
        crash=False
        with open(covfile,'r') as infile:
            for line in infile:
                try:
                    if not line.startswith("#") and line.split("\t")[2] == "mRNA":
                        attrb = line.split("\t")[8]
                        transcript = attrb.split(";")[0].split("=")[1]
                        gene = attrb.split(";")[2].split("=")[1]
                        fully_cov_genes.add(gene)
                        fully_cov_transcripts.add(transcript)
                        self.genesFullyCovered.add(gene)
                        self.transcriptsFullyCovered.add(transcript)
                except:
                    logging.error("Problems parsing the coverage ref file %s" % covfile)
                    logging.info("This step will be skipped.\n")
                    crash=True
                    break

        infile.close()
        if not crash:
            print("Number of genes fully covered by RNA-seq reads in the %s transcriptome.\t%i" % (covfile.split(self.splitstring)[0], len(fully_cov_genes)))
            print("Number of trancripts fully covered by RNA-seq reads in the %s transcriptome.\t%i" % (covfile.split(self.splitstring)[0], len(fully_cov_transcripts)))

    def writeIndividualStats(self,gtf,scaff,genes,newgenes,transcripts,newtranscripts):
        logging.info("Writing individual stats..")
        print("\n####################%s##############" % os.path.basename(gtf).split(self.splitstring)[0])
        print("Total number of scaffolds with genes predicted in the reference\t%i" % len(self.scaffolds))
        print("Number of scaffolds in the transcriptome with expressed transcripts\t%i (%s%%)\n" % (len(scaff), round(len(scaff)/len(self.scaffolds)*100,2)))
        print("Total number of genes in reference\t%i" % len(self.genes))
        print("Total number of genes assembled in this specific transcriptome\t%i" % (len(genes) + len(newgenes)))

        gene_frac_match = round(len(genes)/(len(genes)+len(newgenes))*100,2)
        trans_frac_match = round(len(transcripts)/(len(transcripts)+len(newtranscripts))*100,2)
        self.individualGeneMatchingRate[os.path.basename(gtf).split(self.splitstring)[0]] = gene_frac_match
        self.individualTranscriptMatchingRate[os.path.basename(gtf).split(self.splitstring)[0]] = trans_frac_match

        print("Number of genes in the transcriptome matching the reference\t%i (%s%%)" % (len(genes),gene_frac_match))
        print("Number of unknown new genes\t%i (%s%%)\n" % (len(newgenes),round(len(newgenes)/(len(genes) + len(newgenes))*100,2)))

        print("Total number of transcripts in reference\t%i" % len(self.transcripts))
        print("Total number of transcripts assembled in this specific transcriptome\t%i" % (len(transcripts) + len(newtranscripts)))
        print("Number of transcripts in the transcriptome matching the reference\t%i (%s%%)" % (len(transcripts),trans_frac_match))
        print("Number of unknown new transcripts\t%i (%s%%)\n" % (len(newtranscripts),round(len(newtranscripts)/(len(transcripts) + len(newtranscripts))*100,2)))

    def writeGlobal(self,outprefix):
        with open(outprefix + '_overall-stats.txt', 'w') as outfile:
            outfile.write("Scaffolds in the reference with genes predicted:\t%i\n" % len(self.scaffolds))
            outfile.write("Scaffolds in the set of assembled transcriptomes with expressed transcripts:\t%i\n\n" % len(self.allpresentscaffolds))
            outfile.write("Number of genes in the reference annotation\t%i\n" % len(self.genes))

            overall_fraction_genes = round(len(self.geneInterception)/len(self.genes)*100,2)
            overall_fraction_transcripts = round(len(self.abundance)/len(self.transcripts)*100,2)


            outfile.write("Number of genes from the reference present in the set of assembled transcriptomes\t%i (%s%%)\n\n" % (len(self.geneInterception), overall_fraction_genes))
            outfile.write("Number of transcripts in the reference annotation\t%i\n" % len(self.transcripts))
            outfile.write("Number of transcripts from the reference present in the set of assembled transcriptomes\t%i (%s%%)\n\n" % (len(self.abundance), overall_fraction_transcripts))

            if self.genesFullyCovered:
                outfile.write("Number of genes from the reference fully covered by RNA-seq reads\t%i\n" % len(self.genesFullyCovered))
                outfile.write("Number of transcripts from the reference fully covered by RNA-seq reads\t%i\n\n" % len(self.transcriptsFullyCovered))

            maprate_g,maprate_t = [],[]
            for k,v in self.individualGeneMatchingRate.items():
                maprate_g.append(v)
            for k,v in self.individualTranscriptMatchingRate.items():
                maprate_t.append(v)


            outfile.write("Per sample mean/median mapping rate at gene level\t%s/%s %%\n" % (round(np.mean(maprate_g),2), round(np.median(maprate_g),2)))
            outfile.write("Per sample mean/median mapping rate at transcript level\t%s/%s %%\n\n" % (round(np.mean(maprate_t),2), round(np.median(maprate_t),2)))


            counter = Counter(elem[1] for elem in self.setStrands)
            outfile.write("Forward strand genes matching the reference\t%i\n" % counter["+"])
            outfile.write("Reverse strand genes matching the reference\t%i\n\n" % counter["-"])

        outfile.close()





    def getGlobalStats(self):
        logging.info("Writing global stats..")
        if len(self.individualGeneMatchingRate) == 1:
            logging.info("WARNING.Only one input file was provided, some analysis wont' be performed, and the global stats refer the same as the individual ones.")
            self.writeGlobal(self.outputprefix)
        else:
            self.writeGlobal(self.outputprefix)



def getNumberLines(file):
    return sum(1 for line in open(file))

def main():

    parser = argparse.ArgumentParser(description='Script to validate structural annotation file on the results of transcriptome assembly from RNA-seq data. [e.g. Stringtie].'
                                                 'Starts from the assumption that a good annotation will present very few new assembled transcripts.')
    parser.add_argument(dest='gtf_file', metavar='gtf', nargs="+", help='Annotation file/s to be analyzed. For stringtie assemblies, assumes "ref_gene_id" tag to be present in the atributtes field.')
    parser.add_argument('-s', metavar = 'assembly_software', required= True, help='Software that produced gtf_file. [stringtie]', choices=("stringtie"))
    parser.add_argument('-r', metavar = 'reference_annotation', required = True, help = 'Reference annotation to validate in gff format. Defaults to AUGUSTUS output format.')
    parser.add_argument('-g', metavar = 'genome_table', required = True, help="Genome table file providing information about scaffolds/chromossomes sizes.")
    parser.add_argument('-o', metavar = 'output_basename',required=True,help='Basename to write the output files.')
    parser.add_argument('-l', '--list', action='store_true', help='Input is a file listing all the gtf together, one per line.')
    parser.add_argument('-c', metavar = 'listCoverageRefFiles', help='List of files representing the transcripts from the annotation fully covered by reads. It is an additional'
                                                                    'layer of information produced by stringtie. Number of files must be the same as the gtfs provided and the order must be preserved.')
    parser.add_argument('-i', metavar = 'sampleIdentifier', default= '.' , help ='String to split gtf filenames in order to achieve an unique and shorter sample identifier. Default:"."')
    args = parser.parse_args()

    obj = validateAnnotation(args.s,args.o, args.g, args.i)
    index = 0
    if len(args.gtf_file) > 1 and args.list:
        logging.error("When '-l' set, please don't provide more than 1 file in the positional arguments." )
        exit(1)
    elif args.list and args.c:
        if args.s != "stringtie":
            logging.error("CoverageRef files are only available when stringtie assemblies are produced. Either remove '-c' option or change '-s' to stringtie.")
            exit(1)
        else:
            lines_gtfs = getNumberLines(args.list)
            lines_covRef = getNumberLines(args.c)
            if lines_gtfs != lines_covRef:
                logging.error("Error. Number of gtfs is different than number of coverageRef files. Please check if you provided a coverage ref file in the '-c' argument, instead of a txt file listing those files, one perl line.")
                exit(1)

            obj.processReferenceAnnotGff(args.r)
            with open(args.gtf_file, 'r') as listGtfs:
                with open(args.c, 'r') as covrefs:
                    for line in listGtfs:
                        covfile = ""
                        for l in covrefs:
                            if os.path.abspath(l).rsplit(obj.splitstring, 1)[0] == os.path.abspath(line).rsplit(obj.splitstring,1)[0]:
                                covfile = l.rstrip()
                                break
                        if not covfile:
                            logging.error("Error. File %s does not have a covref file asssociated. Gtfs and CovRef must be the same substring when split by sample identifier ('.', by default)." % line)
                            exit(1)
                        else:
                            obj.processFromStringtie(line, index, lines_gtfs)
                            obj.processCovRefs(covfile)

                    index += 1
                covrefs.close()
            listGtfs.close()
            obj.getGlobalStats()

    elif args.list:
        obj.processReferenceAnnotGff(args.r)
        lines_gtfs = getNumberLines(args.list)
        with open(args.gtf_file, 'r') as listGtfs:
            for line in listGtfs:
                obj.processFromStringtie(line,index, lines_gtfs)
                index += 1
        listGtfs.close()
        obj.getGlobalStats()


    elif args.c:

        obj.processReferenceAnnotGff(args.r)
        numbfiles = len(args.gtf_file)
        lines_covRef = getNumberLines(args.c)
        if numbfiles != lines_covRef:
            logging.error("Error. Number of gtfs is different than number of coverageRef files. Please check if you provided a coverage ref file in the '-c' argument, instead of a txt file listing those files, one perl line.")
            exit(1)

        with open(args.c, 'r') as covrefs:
            for file in args.gtf_file:
                covfile = ""
                for l in covrefs:
                    if os.path.abspath(l).rsplit(obj.splitstring, 1)[0] == os.path.abspath(file).rsplit(obj.splitstring,1)[0]:
                        covfile = l.rstrip()
                        break
                if not covfile:
                    logging.error("Error. File %s does not have a covref file asssociated. Gtfs and CovRef must be the same substring when split by sample identifier ('.', by default)." % file)
                    exit(1)
                elif args.s == "stringtie":
                    obj.processFromStringtie(file, index, numbfiles)
                    obj.processCovRefs(covfile)
                else:
                    logging.error("Ref Coverage files are only available from stringtie assembly.")
                    exit(1)
                index += 1
        covrefs.close()
        obj.getGlobalStats()

    else:
        obj.processReferenceAnnotGff(args.r)
        numbfiles = len(args.gtf_file)
        for file in args.gtf_file:
            if args.s == "stringtie":
                obj.processFromStringtie(file, index, numbfiles)
                index += 1
        obj.getGlobalStats()
if __name__ == "__main__":
    main()