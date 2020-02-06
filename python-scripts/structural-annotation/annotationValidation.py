import matplotlib
matplotlib.use('Agg')
import argparse
import logging
import sys
import os
import numpy as np
from collections import defaultdict
from collections import OrderedDict
from collections import Counter
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
import seaborn as sns
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial import distance
import matplotlib.pyplot as plt



class validateAnnotation():
    def __init__(self, software, outputprefix, genome_table, splitstring, subsampling):
        self.software = software
        self.outputprefix = outputprefix
        self.genome_table = genome_table
        self.splitstring = splitstring
        self.subsampling = subsampling

        self.genomeTableDict = {} #dictionary mapping reference scaffold ids and their length
        self.genesPerScaffold = defaultdict(list) #dictionary listing the existing genes per scaffold in the reference.
        self.genes = set() #set of genes in the reference annotation
        self.transcripts = set() # set of transcripts in the reference annotation
        self.allpresentscaffolds = set() #set of scaffolds with genes predicted in the set of assembled transcripts
        self.abundance = defaultdict(list)  # dictionary of all the transcripts TPM abundance estimation
        self.geneInterception = defaultdict(list) #dictionary of all present/absent genes
        self.individualGeneMatchingRate = OrderedDict() #list representing average of accurate genes predicted
        self.individualTranscriptMatchingRate = OrderedDict()  #list representing average of accurate transcripts predicted
        self.genesFullyCovered = set() #set of genes fully covered by RNA-seq reads (availalbe only when cov ref files are provided)
        self.transcriptsFullyCovered = set() #set of trasncripts fully covered by RNA-seq reads (availalbe only when cov ref files are provided)
        self.setStrands = set() #set of tuples representing the strand of the genes from the reference present in at least on gtf file
        self.indexDict = OrderedDict() #dict mapping file indexes with sample names
        self.genesValidatedPerScaffold = defaultdict(set) #dict listing the validated genes per scaffold

    def processGenomeTable(self,genome):
        logging.info("Processing genome table file..")
        with open(genome, 'r') as infile:
            for line in infile:
                if not line.startswith("#") and not len(line.split("\t")) == 2:
                    logging.error("Please verify genome table file. It should refer to a 2 column tab separated file with contig/scaffold id and its length")
                    exit(1)
                else:
                    id = line.split("\t")[0]
                    length = line.split("\t")[1]
                    self.genomeTableDict[id] = int(length)

    def processReferenceAnnotGff(self,gff):
        logging.info("Processing reference annotation file..")
        with open(gff,'r') as infile:
            for line in infile:
                tmpline = line.rstrip().split("\t")
                if "#" not in line and tmpline[2] == "gene":
                    self.genes.add(tmpline[-1])
                    self.genesPerScaffold[tmpline[0]].append(tmpline[-1])
                elif "#" not in line and tmpline[2] == "transcript":
                    attr = tmpline[-1].split(";")
                    self.transcripts.add(attr[0])

    def processFromStringtie(self,gtf, index, numbfiles):
        logging.info("Reading stringtie transcripts file (%s).." % gtf.rstrip())
        indivscafollds = set()
        indivgenes = set()
        indivnewgenes = set()
        indivtranscript = set()
        indivnewtranscript = set()
        with open(gtf.rstrip(), 'r') as infile:
            self.indexDict[index] = os.path.basename(gtf.rstrip()).split(self.splitstring)[0]
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
                                self.geneInterception[g_id] = [0] * numbfiles
                            self.abundance[t_id][index] = tpm[0].split("\"")[1]
                            self.geneInterception[g_id][index] = 1
                            self.genesValidatedPerScaffold[scaffold_id].add(g_id)

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
            print("Number of genes fully covered by RNA-seq reads in the %s transcriptome.\t%i" % (os.path.basename(covfile.rstrip()).split(self.splitstring)[0], len(fully_cov_genes)))
            print("Number of trancripts fully covered by RNA-seq reads in the %s transcriptome.\t%i\n" % (os.path.basename(covfile.rstrip()).split(self.splitstring)[0], len(fully_cov_transcripts)))

    def writeIndividualStats(self,gtf,scaff,genes,newgenes,transcripts,newtranscripts):
        logging.info("Writing individual stats..")
        print("\n####################%s##############" % os.path.basename(gtf).split(self.splitstring)[0])
        print("Total number of scaffolds with genes predicted in the reference\t%i" % len(self.genesPerScaffold))
        print("Number of scaffolds in the transcriptome with expressed transcripts\t%i\n" % len(scaff))
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
            outfile.write("Scaffolds in the reference with genes predicted:\t%i\n" % len(self.genesPerScaffold))
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
            self.scatterGene2ScaffoldLenght()

            logging.info("Writing set of validated genes and transcripts to file..")
            with open(self.outputprefix + "_validatedGenes.txt", 'w') as outfile:
                for k in self.geneInterception.keys():
                    outfile.write(k + "\n")
            outfile.close()
            with open(self.outputprefix + "_validatedTranscripts.txt", 'w') as outfile:
                for k in self.abundance.keys():
                    outfile.write(k + "\n")
            outfile.close()

        else:
            self.writeGlobal(self.outputprefix)
            self.scatterGene2ScaffoldLenght()
            self.geneMatchingRate()
            self.transcriptMatchingRate()
            self.genePresenceHeatmap()
            self.drawTranscriptAbundanceHeatmap()


            logging.info("Writing whole gene interception set to file..")
            with open(self.outputprefix + "_perSampleGeneComparison.tsv", 'w') as outfile1:
                with open(self.outputprefix + "_validatedGenes.txt", 'w') as outfile2:
                    with open(self.outputprefix + "_validatedTranscripts.txt", 'w') as outfile3:
                        outfile1.write("#gene_id" + "\t" + '\t'.join(self.indexDict.values()) + "\n")
                        for k,v in self.geneInterception.items():
                            outfile1.write(k + "\t" + '\t'.join(str(i) for i in v) + "\n")
                            outfile2.write(k + "\n")
                        for k in self.abundance.keys():
                            outfile3.write(k + "\n")
            outfile1.close()
            outfile2.close()
            outfile3.close()

            logging.info("Writing transcripts abundance interception to file..")
            with open(self.outputprefix + "_transcriptAbundance.tsv", 'w') as outfile:
                outfile.write("#transcript_id" + "\t" + '\t'.join(self.indexDict.values()) + "\n")
                for k,v in self.abundance.items():
                    outfile.write(k + "\t" + '\t'.join(str(i) for i in v) + "\n")
            outfile.close()

    def scatterGene2ScaffoldLenght(self):

        logging.info("Drawing scatter plot about gene validation in function of scaffold length..")
        fig,ax = plt.subplots()


        present_scaf_length, absent_scaf_length, n_present, n_absent = [],[],[],[]
        genesNotValidated = defaultdict(set)
        for scaff_id, scaff_length in sorted(self.genomeTableDict.items()):
            if scaff_id in self.genesPerScaffold:

                if scaff_id in self.genesValidatedPerScaffold.keys():
                    present_scaf_length.append(scaff_length)
                    n_present.append(len(self.genesValidatedPerScaffold[scaff_id]))

                for gene_id in self.genesPerScaffold[scaff_id]:

                    if not scaff_id in self.genesValidatedPerScaffold.keys() or not gene_id in self.genesValidatedPerScaffold[scaff_id]:
                        genesNotValidated[scaff_id].add(gene_id)

        for k,v in genesNotValidated.items():
            absent_scaf_length.append(self.genomeTableDict[k])
            n_absent.append(len(v))


        logging.info("Scaffolds with genes predicted\t%i "% len(self.genesPerScaffold))
        logging.info("Scaffolds with genes predicted and validated\t%i" % len(self.genesValidatedPerScaffold))
        logging.info("Scaffolds with genes predicted and not validated by RNA-seq\t%i" % len(genesNotValidated))
        logging.info("Validated genes\t%i" % sum(len(v) for v in self.genesValidatedPerScaffold.values()))
        logging.info("Invalidated genes\t%i" % sum(len(v) for v in genesNotValidated.values()))


        fig, ax1 = plt.subplots()
        ax1.scatter(present_scaf_length, n_present, color ='olivedrab')
        ax1.set_xlabel('Scaffolds length (bp)')
        ax1.set_ylabel('Genes predicted and validated', color='olivedrab')
        for tl in ax1.get_yticklabels():
            tl.set_color('olivedrab')

        ax2 = ax1.twinx()
        ax2.scatter(absent_scaf_length, n_absent, color = 'maroon')
        ax2.set_ylabel('Genes predicted and not validated', color='maroon')
        for tl in ax2.get_yticklabels():
            tl.set_color('maroon')

        ax1.set_xlim(min(self.genomeTableDict.values()) ,max(self.genomeTableDict.values()) + 10000)
        ax1.set_ylim(0,max(n_present) + 10)
        ax2.set_xlim(min(self.genomeTableDict.values()) ,max(self.genomeTableDict.values()) + 10000)
        ax2.set_ylim(0,max(n_absent) + 10)
        plt.savefig(self.outputprefix + "_geneValidation_allScaffolds.png")
        plt.close()

        logging.info("Drawing subplot..")
        fig2,ax3=plt.subplots()
        ax3.scatter(present_scaf_length, n_present, color ='olivedrab')
        ax3.set_xlabel('Scaffolds length (bp)')
        ax3.set_ylabel('Genes predicted and validated', color='olivedrab')
        for tl in ax3.get_yticklabels():
            tl.set_color('olivedrab')

        ax4 = ax3.twinx()
        ax4.scatter(absent_scaf_length, n_absent, color = 'maroon')
        ax4.set_ylabel('Genes predicted and not validated', color='maroon')
        for tl in ax4.get_yticklabels():
            tl.set_color('maroon')

        ax3.set_xlim(1000,20000)
        ax3.set_ylim(-1,20)
        ax4.set_xlim(1000,20000)
        ax4.set_ylim(-1,20)
        plt.savefig(self.outputprefix + "_geneValidation_smallScaffolds.png")
        plt.close()

    def genePresenceHeatmap(self):
        Cols=list(self.indexDict.values())

        plt.figure()
        logging.info("Drawing gene expression presence/absence heatpmap..")
        df = pd.DataFrame.from_dict(self.geneInterception, orient='index')
        df.columns=Cols
        df = df[df.columns].astype(int)

        logging.info("\tExpressed genes in reference %i" % len(df))
        #Remove row with all samples having the gene expressed
        df = df[(df == 0).any(axis=1)]
        logging.info("\tGenes with at least one sample with no expression %i" % len(df))
        #sort and subsample first N rows
        #df = df.sort_values(by=Cols,axis=0, ascending=False)[:self.subsampling]

        #random subsample of N rows
        if df.shape[0] > self.subsampling:
            df = df.sample(n=self.subsampling)
        #df2=df[:200]
        else:
            logging.info("Subsampling provided is larger than all the rows in dataframe. Whole df will be used (%i)" % df.shape[0])

        #row_linkage = hierarchy.linkage(distance.pdist(df), method='average')
        logging.info("\tCreating col linkage")
        col_linkage = hierarchy.linkage(distance.pdist(df.T), method='average')
        p = sns.clustermap(df,col_linkage=col_linkage,figsize=(20,25), linewidths=0, yticklabels=False, metric="correlation")


        #plt.setp(p.ax_heatmap.get_yticklabels(), rotation=0)
        plt.savefig(self.outputprefix + "_genesPresence.png")
        plt.close()

    def drawTranscriptAbundanceHeatmap(self):
        Cols = list(self.indexDict.values())

        plt.figure()
        logging.info("Drawing transcript abundances heatmap..")

        df = pd.DataFrame.from_dict(self.abundance, orient='index')
        df.columns = Cols
        df = df[df.columns].astype(float)
        df = df.fillna(0)

        logging.info("\tExpressed transcripts in reference %i" % len(df))
        #select row with at least 1 sample with 0 expression
        df = df[(df == 0).any(axis=1)]
        logging.info("\tTranscripts with at least on sample with no expression %i" % len(df))
        #df = df.nlargest(self.subsampling,df.sum(axis=1))

        if df.shape[0] > self.subsampling:
            df['sum'] = df.sum(axis=1)
            df = df.nlargest(self.subsampling,'sum').drop('sum',1)
        #df2=df[:200]
        else:
            logging.info("Subsampling provided is larger than all the rows in dataframe. Whole df will be used (%i)" % df.shape[0])


        logging.info("\tCreating col linkage ")
        col_linkage = hierarchy.linkage(distance.pdist(df.T), method='average')
        #logging.info("Createing row linkage")
        #row_linkage = hierarchy.linkage(distance.pdist(df), method='average')


        cmap = sns.cubehelix_palette(as_cmap=True, rot=-.3, light=1)
        p = sns.clustermap(df,col_linkage=col_linkage,figsize=(20,25), standard_scale=0,linewidths=0, yticklabels=False,metric="correlation", cmap=cmap)
        #p = sns.heatmap(df,cmap=cmap)

        #plt.tight_layout()
        plt.setp(p.ax_heatmap.get_yticklabels(), visible=False,rotation=0)
        plt.yticks(visible=False)
        plt.title("Top %i transcripts abundance." % self.subsampling, loc='center',fontsize=8)
        plt.savefig(self.outputprefix + "_transcriptsAbundance.png",bbox_inches='tight',dpi=100)
        plt.close()


    def geneMatchingRate(self):
        logging.info("Drawing individual sample gene matching rate..")
        plt.figure(figsize=(12,9))
        plt.bar(range(len(self.individualGeneMatchingRate)), self.individualGeneMatchingRate.values())
        plt.xticks(range(len(self.individualGeneMatchingRate)), self.individualGeneMatchingRate.keys(), rotation='vertical')
        plt.title("Per sample fraction of expressed genes that match the reference annotation")
        plt.xlabel("Sample")
        plt.ylabel("Fraction (%)")
        plt.rcParams.update({'font.size': 12})
        plt.tight_layout()
        plt.savefig(self.outputprefix + "_geneBarplot.png")
        plt.close()

    def transcriptMatchingRate(self):
        logging.info("Drawing individual sample transcript matching rate..")
        plt.figure(figsize=(12,9))
        plt.bar(range(len(self.individualTranscriptMatchingRate)), self.individualTranscriptMatchingRate.values())
        plt.xticks(range(len(self.individualTranscriptMatchingRate)), self.individualTranscriptMatchingRate.keys(), rotation='vertical')
        plt.title("Per sample fraction of expressed transcripts that match the reference annotation")
        plt.xlabel("Sample")
        plt.ylabel("Fraction (%)")
        plt.rcParams.update({'font.size': 12})
        plt.tight_layout()
        plt.savefig(self.outputprefix + "_transcriptBarplot.png")
        plt.close()

def main():

    parser = argparse.ArgumentParser(description='Script to validate structural annotation file on the results of transcriptome assembly from RNA-seq data. [e.g. Stringtie].'
                                                 'Starts from the assumption that a good annotation will present very few new assembled transcripts.')
    parser.add_argument(dest='gtf_file', metavar='gtf', nargs="+", help='Annotation file/s to be analyzed. For stringtie assemblies, assumes "ref_gene_id" tag to be present in the atributtes field.')
    parser.add_argument('-p', metavar = 'assembly_software', required= True, help='Software that produced gtf_file. [stringtie]', choices=("stringtie"))
    parser.add_argument('-r', metavar = 'reference_annotation', required = True, help = 'Reference annotation to validate in gff format. Defaults to AUGUSTUS output format.')
    parser.add_argument('-g', metavar = 'genome_table', required = True, help="Genome table file providing information about scaffolds/chromossomes sizes.")
    parser.add_argument('-o', metavar = 'output_basename',required=True,help='Basename to write the output files.')
    parser.add_argument('-l', '--list', action='store_true', help='Input is a file listing all the gtf together, one per line.')
    parser.add_argument('-c', metavar = 'listCoverageRefFiles', help='List of files representing the transcripts from the annotation fully covered by reads. It is an additional'
                                                                    'layer of information produced by stringtie. Number of files must be the same as the gtfs provided and the order must be preserved.')
    parser.add_argument('-i', metavar = 'sampleIdentifier', default= '.' , help ='String to split gtf filenames in order to achieve an unique and shorter sample identifier. Default:"."')
    parser.add_argument('-s', metavar = 'subsampling', default=1000, type=int, help ='Number of features (transcripts and genes) to be subsampled when drawing clustered heatmaps. This is only useful when more than 1 gtf file is provided.'
                                                                                     ' Please be carefull that large values will make the whole process extremely slow.Default:1000')
    args = parser.parse_args()

    obj = validateAnnotation(args.p,args.o, args.g, args.i, args.s)
    index = 0
    if len(args.gtf_file) > 1 and args.list:
        logging.error("When '-l' set, please don't provide more than 1 file in the positional arguments." )
        exit(1)
    elif args.list and args.c:

        if args.p != "stringtie":
            logging.error("CoverageRef files are only available when stringtie assemblies are produced. Either remove '-c' option or change '-s' to stringtie.")
            exit(1)
        else:
            lines_gtfs = sum(1 for line in open(args.gtf_file[0]))
            lines_covRef = sum(1 for line in open(args.c))
            if lines_gtfs != lines_covRef:
                logging.error("Error. Number of gtfs is different than number of coverageRef files. Please check if you provided a coverage ref file in the '-c' argument, instead of a txt file listing those files, one perl line.")
                exit(1)
            obj.processGenomeTable(args.g)
            obj.processReferenceAnnotGff(args.r)
            with open(args.gtf_file[0], 'r') as listGtfs:
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
        obj.processGenomeTable(args.g)
        obj.processReferenceAnnotGff(args.r)
        lines_gtfs = sum(1 for line in open(args.gtf_file[0]))
        if args.p == "stringtie":
            with open(args.gtf_file[0], 'r') as listGtfs:
                for line in listGtfs:
                    obj.processFromStringtie(line,index, lines_gtfs)
                    index += 1
            listGtfs.close()
            obj.getGlobalStats()


    elif args.c:
        obj.processGenomeTable(args.g)
        obj.processReferenceAnnotGff(args.r)
        numbfiles = len(args.gtf_file)
        lines_covRef = sum(1 for line in open(args.c))
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
                elif args.p == "stringtie":
                    obj.processFromStringtie(file, index, numbfiles)
                    obj.processCovRefs(covfile)
                else:
                    logging.error("Ref Coverage files are only available from stringtie assembly.")
                    exit(1)
                index += 1
        covrefs.close()
        obj.getGlobalStats()

    else:
        obj.processGenomeTable(args.g)
        obj.processReferenceAnnotGff(args.r)
        numbfiles = len(args.gtf_file)
        for file in args.gtf_file:
            if args.p == "stringtie":
                obj.processFromStringtie(file, index, numbfiles)
                index += 1
        obj.getGlobalStats()
if __name__ == "__main__":
    main()
