import argparse
import logging
import os
import sys
from collections import defaultdict
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
from tempfile import NamedTemporaryFile
import codecs
from pybedtools import BedTool
from operator import itemgetter
import numpy as np
def processGenomeTable(genome):
    genomeTableDict={}
    logging.info("Processing genome table file..")
    with open(genome, 'r') as infile:
        for line in infile:
            if not line.startswith("#") and not len(line.split("\t")) == 2:
                logging.error("Please verify genome table file. It should refer to a 2 column tab separated file with contig/scaffold id and its length")
                exit(1)
            else:
                id = line.split("\t")[0]
                length = line.split("\t")[1]
                genomeTableDict[id] = int(length)
    return  genomeTableDict

def processM4(file):

    queryDict,refDict=defaultdict(list),defaultdict(list)
    with open(file,'r') as m4file:
        for line in m4file:
            if not line.startswith("#") or not line.startswith("qName"):
                (query,refcontig,score,identityStr,qstrand,qalignstartStr,qalignendStr,qlengthStr,rstrand,ralignstartStr,ralignendStr,rlengthStr,mapQStr) = line.rstrip().split()
                identity=float(identityStr)
                qalignstart=int(qalignstartStr)
                qalignend=int(qalignendStr)
                qlength=int(qlengthStr)
                ralignstart=int(ralignstartStr)
                ralignend=int(ralignendStr)
                rlength=int(rlengthStr)
                mapQ=int(mapQStr)
                if "|" in refcontig:
                    ref_contig = refcontig.split("|")[0]
                else:
                    ref_contig=refcontig	
                if qstrand == 1: #reverse, qstrand field will be kept the same to distinguish from original reverse alignments
                    qnewstart=qlength-qalignend
                    qneweend=qlength-qnewstart
                    queryDict[query].append((ref_contig,score, identity, qstrand,qnewstart,qneweend,qlength, mapQ))
                else:
                    queryDict[query].append((ref_contig,score, identity, qstrand,qalignstart,qalignend,qlength, mapQ))


                if rstrand == 1:
                    rnewstart=rlength-ralignend
                    rnewend=rlength-rnewstart
                    refDict[ref_contig].append((identity,rstrand,rnewstart,rnewend, mapQ))
                else:
                    refDict[ref_contig].append((identity,rstrand,ralignstart,ralignend, mapQ))

    return queryDict,refDict

def getStats(queryDict,refDict,refLength,basename,noCoverage,genomeTblFile):

    alignReadFraction, identityQuery, mapQ=[],[],[]
    multipleMappedReads,multipleReferenceMaps,fakeMultipleMap=defaultdict(list),defaultdict(set),[]
    identityPerScaff=defaultdict(list)
    logging.info("Processing queries..")
    for query,align in queryDict.items():
        for al in align:
            alnFraction=round((al[5]-al[4])/al[6]*100,2)
            if alnFraction > 100:
                print(query, align)
            if len(align) > 1:#multiple alignment

                uniqReadLen=set([e[6] for e in align]) #get set of unique sequences lengths
                if len(uniqReadLen) > 1:
                    fakeMultipleMap.append(query)
                    logging.warning("%s is a duplicate ID. Original read ID may have white spaces, but as blasr trimmed read name up to the first one found, they acount for the same. Read discarded to acount for the mapping stats." % query)
                else:
                    multipleReferenceMaps[query].add(al[0])
                    multipleMappedReads[query].append((al[2],alnFraction))

            alignReadFraction.append((alnFraction))
            identityQuery.append(al[2])
            mapQ.append(al[7])


    logging.info("Processing reference..")
    for contig, align in refDict.items():
        for al in align:
            identityPerScaff[contig].append(al[0])

    with open(basename + "_stats.txt", 'w') as outfile:
        outfile.write("Number of reads aligned:\t%i\n" % len(queryDict))
        outfile.write("Number of reads with one unique alignment [only makes sense if --bestn argument was set to be != 1]:\t%i\n" % (len(queryDict) - len(multipleMappedReads) - len(fakeMultipleMap)))
        outfile.write("Number of true multiple alignments [includes split read alignments, makes sense in --bestn != 1]:\t%i\n" % len(multipleMappedReads))
        outfile.write("Number of problematic alignment because of duplicates:\t%i\n\n" % len(fakeMultipleMap))
        outfile.write("Average read mapping identity:\t%f\n" % np.mean(identityQuery))
        outfile.write("Minimum/Maximum read mapping identity:\t%f/%f\n" % (min(identityQuery),max(identityQuery)))
        outfile.write("Per alignment average fraction of the reads length aligned:\t%f\n" % np.mean(alignReadFraction))
        outfile.write("Average mapping quality:\t%f\n" % round(np.mean(mapQ),2))

        #improve code here
        multipleContigs=[]
        for query, contigs in multipleReferenceMaps.items():
            if len(contigs) > 1:
                multipleContigs.append(query)
        outfile.write("Number of reads belonging to the multiple alignments class mapping to different contigs in the referece:\t%i\n\n" % len(multipleContigs))

        outfile.write("Number of scaffolds in reference:\t%i\n" % len(refLength))
        outfile.write("Number of reference scaffolds with alignments:\t%i\n" % len(refDict))
        outfile.close()
        if noCoverage:
            logging.info("Starting genome coverage analysis..")
            iterable = createTmpBedFiles(refDict,genomeTblFile)
            getGenomeCov(iterable,refDict,refLength,basename)

        with open(basename + "_stats.txt", 'a') as outfile:
            outfile.write("\n#Per scaffold average identity of read alignments:\n")
            for contig, identities in identityPerScaff.items():
                outfile.write("%s\t%f\n" % (contig,round(np.mean(identities),2)))


def createTmpBedFiles(refDict, refGenomeTable):
    logging.info("Creating tmp bed file..")
    tmpfileRef = NamedTemporaryFile(delete=True)
    i=1
    with codecs.open(tmpfileRef.name, 'w+b', encoding='utf-8') as tmp1:
            for key,val in refDict.items():
                for tuple in sorted(val, key=itemgetter(2)):
                    name='l'+str(i)
                    tmp1.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (key,tuple[2],tuple[3],name,0,'+'))
                    i+=1
    tmp1.close()
    a=BedTool(tmpfileRef.name)
    logging.info("Merging sets of intervals and calculating genome coverage for reference genome..")
    iterable_ref = a.merge().genome_coverage(g=refGenomeTable)
    return iterable_ref

def getGenomeCov(iterable_ref,refDict,refGeneTable,outbasename):
    finalREF=defaultdict(list)
    unl_ref=[]

    with open(outbasename + "_stats.txt", 'a') as outfile:
        outfile.write("\n#Genome coverage analysis#\n")

        for line in iterable_ref:
            attr=str(line).rstrip().split("\t")

            if attr[0] != 'genome' and attr[1] == "1":
                finalREF[attr[0]] = [attr[3],attr[2],attr[4]]
            elif attr[0] == 'genome' and attr[1] == "1":
                outfile.write("%s\t%i\n" % ("Total length of reference genome:",int(attr[3])))
                outfile.write("%s\t%i%s%f%s\n" % ("Fracion of reference genome aligned:",int(attr[2])," [",float(attr[4]),"]"))
                unl_len=0
                for k in refGeneTable.keys():
                    if not k in refDict.keys():
                        unl_len+=int(refGeneTable[k])
                        unl_ref.append((k,refGeneTable[k]))
                outfile.write("%s\t%i%s%f%s\n" % ("Total length considering the fully unaligned scaffolds:",unl_len," [",round(unl_len/int(attr[3]),4), "]"))
    outfile.close()
    with open(outbasename + "_referenceGenomeCov.tsv",'w') as outfileref:
        outfileref.write("#%s\t%s\t%s\t%s\n" % ('scaffold_id','scaffold_length','aligned_bases','fraction_covered'))
        for k,v in finalREF.items():
            outfileref.write("%s\t%s\n" % (k,'\t'.join(v)))

        for k,v in sorted(unl_ref, reverse=True):
            outfileref.write("%s\t%i\t%i\t%i\n" % (k,v,0,0))
    outfileref.close()


def main():

    parser = argparse.ArgumentParser(description='Script to analyse overall stats from an alignment of pacbio reads to a reference.')
    parser.add_argument(dest='input_file', metavar='blasrFile', nargs="+", help='Alignment file/s to be processed.')
    parser.add_argument('-g', metavar = 'genome_table', required = True, help="Genome table file providing information about reference scaffolds/chromossomes sizes.")
    parser.add_argument('-l', '--list', action='store_true', help='Input is a file listing all the alignment files together, one per line.')
    parser.add_argument('-f','--format', default='m4',help='Format of blasr output files. Default:[m4].')
    parser.add_argument('-n','--noCoverage', action='store_false',help='Flag to disable reference coverage analysis. Default: enabled!')
    args = parser.parse_args()

    if len(args.input_file) > 1 and args.list:
        logging.info("When '-l' set, please don't provide more than 1 file in the positional arguments." )
        exit(1)

    elif args.list:
        with open(args.input_file[0], 'r') as listBlasrOut:
            for line in listBlasrOut:
                l=line.rstrip()
                basename = os.path.splitext(os.path.basename(l))[0]
                refGeneTable=processGenomeTable(args.g)
                logging.info("Processing %s alignment file.." % (basename))
                queryDict,refDict=processM4(l)
                getStats(queryDict,refDict,refGeneTable,basename,args.noCoverage,args.g)


    else:
        for file in args.input_file:
            basename=os.path.splitext(os.path.basename(file))[0]
            refGeneTable=processGenomeTable(args.g)
            logging.info("Processing %s alignment file.." % (basename))
            queryDict,refDict=processM4(file)
            getStats(queryDict,refDict,refGeneTable,basename,args.noCoverage,args.g)


if __name__ == "__main__":
    main()
