import argparse
from collections import OrderedDict
from collections import defaultdict
from tempfile import NamedTemporaryFile
import codecs
from pybedtools import BedTool
from collections import Counter

def createGenomeTableDict(refGenomeTable,queryGenomeTable):
    refLenghtsDict,queryLengthDict=OrderedDict(), OrderedDict()
    with open(refGenomeTable, 'r') as infile1:
        for line in infile1:
            if len(line.split("\t")) != 2:
                print("Error, genome table files must have exactly 2 columns.")
                exit(1)
            elif line.split("\t")[0] in refLenghtsDict:
                print("Error, repeated scaffold id in genome table file %s" % refGenomeTable)

            else:
                refLenghtsDict[line.split("\t")[0]] = line.split("\t")[1].rstrip()
    infile1.close()


    with open(queryGenomeTable, 'r') as infile2:
        for line in infile2:
            if len(line.split("\t")) != 2:
                print("Error, genome table files must have exactly 2 columns.")
                exit(1)
            elif line.split("\t")[0] in queryLengthDict:
                print("Error, repeated scaffold id in genome table file %s" % queryGenomeTable)

            else:
                queryLengthDict[line.split("\t")[0]] = line.split("\t")[1].rstrip()
    infile2.close()
    return refLenghtsDict,queryLengthDict

def processLastal(lastOut, isAllvsAll, removeIsoformAln):
    alignedRef,alignedQuery, idsMatchDict = defaultdict(list), defaultdict(list), defaultdict(set)
    with open(lastOut,"r") as infile:

        for l in infile:
            line=l.rstrip()
            if not line.startswith("#") and len(line.split()) > 10:
                (score, scaff_reference, refstart, refAlignSize, refstrand, scaffSize, scaff_query, qstart, qAlignSize, qstrand, qScaffSize, blocks, mismatch) = line.split()[:13]

   #             if isAllvsAll and removeIsoformAln and scaff_reference.rsplit('.',1)[0] == scaff_query.rsplit('.',1)[0]:
   #                 allType+=1

                if isAllvsAll and removeIsoformAln and scaff_reference.rsplit('.',1)[0] != scaff_query.rsplit('.',1)[0]:
                    alignedRef[scaff_reference].append((int(refstart),int(refstart)+int(refAlignSize)-1))
                    alignedQuery[scaff_query].append((int(qstart),int(qstart)+int(qAlignSize)-1,qstrand))
                    idsMatchDict[scaff_reference].add(scaff_query)
   #             elif isAllvsAll  and scaff_reference != scaff_query and scaff_reference.rsplit('.',1)[0] == scaff_query.rsplit('.',1)[0]:
   #                 isoformFault+=1

   #             elif isAllvsAll and scaff_reference == scaff_query:
   #                 selfMatch+=1

                elif isAllvsAll and not removeIsoformAln and scaff_reference != scaff_query:
                    alignedRef[scaff_reference].append((int(refstart),int(refstart)+int(refAlignSize)-1))
                    alignedQuery[scaff_query].append((int(qstart),int(qstart)+int(qAlignSize)-1,qstrand))
                    idsMatchDict[scaff_reference].add(scaff_query)

                elif not isAllvsAll:
                    alignedRef[scaff_reference].append((int(refstart),int(refstart)+int(refAlignSize)-1))
                    alignedQuery[scaff_query].append((int(qstart),int(qstart)+int(qAlignSize)-1,qstrand))

    return alignedRef,alignedQuery,idsMatchDict

def writeMatchesIDs(idsMatch,basename):
    with open(basename + "_IDsAligned.txt", 'w') as out:
        for k,v in idsMatch.items():
            out.write('%s\t%s\t%i\n' % (k, ','.join(str(s) for s in v), len(v)))

def transformNegativeStrandCoordinates(alignedQuery,dictQueryGenomeTable):
    #Lastal reports matches in the negative strand of the query by displaying reverse completements coordinates of it.
    #Therefore, for all of the query scaffolds for which we have a mix on the alignments strands, we will transform the negative ones in their
    #reverse complement of the positive strand, so we can access true alignment coverage in such cases.
    #In scaffolds for which the aligned blocks represent negative strand alignents, we will mantain the configuration, as we are just interested
    #in measuring the ammount of the query genome that aligned.

    transformedQuerydict=defaultdict(list)
    mixedAlnFraction= defaultdict(list)
    totallyAlnForwd,totallyAlnRev=[],[]

    for k,v in alignedQuery.items():

        number_blocks = len(v)
        counter = Counter(elem[2] for elem in v)
        if len(counter) == 1 and number_blocks == counter['+']:
            totallyAlnForwd.append((k,counter['+']))
            transformedQuerydict[k] = v
        elif len(counter) == 1 and number_blocks == counter['-']:
            totallyAlnRev.append((k,counter['-']))
            transformedQuerydict[k] = v
        else:
            pos=round(counter['+']/number_blocks,4)
            neg=round(counter['-']/number_blocks,4)
            mixedAlnFraction[k] = [counter['+'],counter['-'],pos,neg]

            finalCoords=[]
            for block in v:
                if block[2] == '-':
                    endCoord=int(dictQueryGenomeTable[k]) - block[0] - 1
                    beginCoord=int(dictQueryGenomeTable[k]) - block[1] - 1
                    finalCoords.append((beginCoord,endCoord,'+','neg'))
                else:
                    finalCoords.append((block[0],block[1],'+'))

            transformedQuerydict[k] = sorted(finalCoords, key=lambda x: x[0])

    return transformedQuerydict,totallyAlnForwd,totallyAlnRev, mixedAlnFraction


def createtmpBedFiles(alnRef,alntasnformedQuery,referenceGenomeTable,alignedGenomeTable):
    tmpfileRef = NamedTemporaryFile(delete=True)
    tmpfileQuery = NamedTemporaryFile(delete=True)
    i,j=1,1
    with codecs.open(tmpfileRef.name, 'w+b', encoding='utf-8') as tmp1:
        with codecs.open(tmpfileQuery.name, 'w+b', encoding='utf-8') as tmp2:
            for key,val in alnRef.items():
                for tuple in sorted(val, key=lambda x: x[0]):
                    name='line'+str(i)
                    tmp1.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (key,tuple[0],tuple[1],name,0,'+'))
                    i+=1

            for key,val in alntasnformedQuery.items():
                for tuple in sorted(val, key=lambda x: x[0]):
                    name='line'+str(j)
                    tmp2.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (key,tuple[0],tuple[1],name,0,'+'))
                    j+=1
    tmp1.close()
    tmp2.close()
    a=BedTool(tmpfileRef.name)
    b=BedTool(tmpfileQuery.name)
    print("Merging sets of intervals and calculating genome coverage for reference genome..")
    iterable_ref = a.merge().genome_coverage(g=referenceGenomeTable)
    print("Merging sets of intervals and calculating genome coverage for aligned genome..")
    iterable_query = b.merge().genome_coverage(g=alignedGenomeTable)

    return iterable_ref,iterable_query

def generateStats(it_ref,it_query,refDict, queryDict,refGenTable,querGenTable,forwdAlnScaf,revAlnScaf,mixAlnScaf, outBasename):

    finalREF=defaultdict(list)
    finalQUERY=defaultdict(list)
    unl_query,unl_ref=[],[]
    with open(outBasename + "_stats.txt", 'w') as outstats:
        for line in it_ref:
            attr=str(line).rstrip().split("\t")

            if attr[0] != 'genome' and attr[1] == "1":
                finalREF[attr[0]] = [attr[3],attr[2],attr[4],str(len(refDict[attr[0]]))]
            elif attr[0] == 'genome' and attr[1] == "1":
                outstats.write("%s\t%i\n" % ("Total length of reference genome:",int(attr[3])))
                outstats.write("%s\t%i%s%f%s\n" % ("Fracion of reference genome aligned:",int(attr[2])," [",float(attr[4]),"]"))
                outstats.write("%s\t%i\n" % ("Number of scaffolds in reference:", len(refGenTable)))
                outstats.write("%s\t%i\n" % ("Number of reference scaffolds with alignments:", len(refDict)))
                unl_len=0
                for k in refGenTable.keys():
                    if not k in refDict.keys():
                        unl_len+=int(refGenTable[k])
                        unl_ref.append((k,refGenTable[k]))
                outstats.write("%s\t%i%s%f%s\n\n\n" % ("Total length considering the fully unaligned scaffolds:",unl_len," [",round(unl_len/int(attr[3]),4), "]"))

        for line in it_query:
            attr=str(line).rstrip().split("\t")
            if attr[0] != 'genome' and attr[1] == "1":
                finalQUERY[attr[0]] = [attr[3],attr[2],attr[4],str(len(queryDict[attr[0]]))]
            elif attr[0] == 'genome' and attr[1] == "1":
                outstats.write("%s\t%i\n" % ("Total length of query genome:",int(attr[3])))
                outstats.write("%s\t%i%s%f%s\n" % ("Fracion of query genome aligned:",int(attr[2])," [",float(attr[4]),"]"))
                outstats.write("%s\t%i\n" % ("Number of scaffolds in query:", len(querGenTable)))
                outstats.write("%s\t%i\n" % ("Number of query scaffolds with alignments:", len(queryDict)))
                unl_len=0
                for k in querGenTable.keys():
                    if not k in queryDict:
                        unl_len+=int(querGenTable[k])
                        unl_query.append((k,querGenTable[k]))
                outstats.write("%s\t%i%s%f%s\n" % ("Total length considering the fully unaligned scaffolds:",unl_len," [",round(unl_len/int(attr[3]),4), "]"))
                outstats.write("%s\t%i\n" % ("Number of query scaffolds with all aligned blocks in the forward strand:", len(forwdAlnScaf)))
                outstats.write("%s\t%i\n" % ("Number of query scaffolds with all aligned blocks in the reverse strand:", len(revAlnScaf)))
                outstats.write("%s\t%i\n" % ("Number of query scaffolds with aligned blocks in both strands:", len(mixAlnScaf)))

    outstats.close()
    for tup in forwdAlnScaf:
        finalQUERY[tup[0]].append(str(tup[1]))
        finalQUERY[tup[0]].append("0")
    for tup in revAlnScaf:
        finalQUERY[tup[0]].append("0")
        finalQUERY[tup[0]].append(str(tup[1]))
    for k, v in mixAlnScaf.items():
        finalQUERY[k].append(str(v[0]))
        finalQUERY[k].append(str(v[1]))

    with open(outBasename + "_referenceGenome.tsv",'w') as outfileref:
        outfileref.write("#%s\t%s\t%s\t%s\t%s\n" % ('scaffold_id','scaffold_length','aligned_bases','fraction_covered','last_blocksAligned'))
        for k,v in finalREF.items():
            outfileref.write("%s\t%s\n" % (k,'\t'.join(v)))
    outfileref.close()

    with open(outBasename + "_queryGenome.tsv", 'w') as outfilequery:
        outfilequery.write("#%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ('scaffold_id','scaffold_length','aligned_bases','fraction_covered','last_blocksAligned',
                                                              'forwardAligned_blocks','reverseAligned_blocks'))
        for k,v in finalQUERY.items():
            outfilequery.write("%s\t%s\n" % (k,'\t'.join(v)))
    outfilequery.close()

    with open(outBasename + "_unalignedScaffolds_reference.txt", 'w') as out:
        for elem in unl_ref:
            out.write(elem[0] + "\t" + elem[1] + "\n")
    out.close()

    with open(outBasename + "_unalignedScaffolds_query.txt", 'w') as out:
        for elem in unl_query:
            out.write(elem[0] + "\t" + elem[1] + "\n")
    out.close()

def main():
    parser = argparse.ArgumentParser(description='Script to analyse lastal tab output file.')
    parser.add_argument(dest='lastOutput', metavar='lastTabOut', help='Lastal output file in TAB format mapping the reference inputted to the query.')
    parser.add_argument(dest='referenceGenomeTable', metavar='refGenomeTable', help='Tab delimited file displaying reference sequences and their length.')
    parser.add_argument(dest='alignedGenomeTable', metavar='alnGenomeTable',  help='Tab delimited file dsiplaying query sequences and their length.')
    parser.add_argument(dest='outputBasemame', metavar='outputBasename', help='Basename to write output files.')
    parser.add_argument('-a','--all2all', action = 'store_true', help='Lastal run refers to an all vs all alignment of the reference sequences. If so, self sequence '
                                                                      'matches will not be processed.')
    parser.add_argument('-i','--isoformsRemoval', action = 'store_true', help='When comparing transcript sequences of one single reference [-a argument], this option '
                                                                              'will also not process matches of different isoforms of the same gene. Script detects a different '
                                                                              'isoform when reference and query names share the same string up to the last \'.\'. E.g gene1.t1;'
                                                                              'gene1.t2 will be treated as a different isoform for the same gene.')
    args = parser.parse_args()

    if args.isoformsRemoval and args.all2all is None:
        parser.error("-i argument requires '-a'.")
        exit(1)

    dictRef,dictQuery=createGenomeTableDict(args.referenceGenomeTable, args.alignedGenomeTable)
    print("Processing lastal output..")
    alignRef,alignQuery,idsMatchDict=processLastal(args.lastOutput,args.all2all,args.isoformsRemoval)
    print("Reporting alignments IDs")
    writeMatchesIDs(idsMatchDict, args.outputBasemame)
    print("Transforming negative strand alignments")
    transformedQueryDic,totallyAlnForwd,totallyAlnRev,mixedAlnFraction= transformNegativeStrandCoordinates(alignQuery,dictQuery)
    print("Creating tmp bed files ..")
    it_ref,it_query = createtmpBedFiles(alignRef,transformedQueryDic, args.referenceGenomeTable, args.alignedGenomeTable)
    print("Generating stats..")
    generateStats(it_ref,it_query,alignRef,transformedQueryDic,dictRef,dictQuery,totallyAlnForwd,totallyAlnRev,mixedAlnFraction,args.outputBasemame)
    print("Done!")


if __name__ == "__main__":
    main()
