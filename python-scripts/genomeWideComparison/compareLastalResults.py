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

def processLastal(lastOut, dictRef, dictQuery):
    alignedRef,alignedQuery = defaultdict(list), defaultdict(list)
    with open(lastOut,"r") as infile:

        for l in infile:
            line=l.rstrip()
            if not line.startswith("#") and len(line.split()) > 10:
                (score, scaff_reference, refstart, refAlignSize, refstrand, scaffSize, scaff_query, qstart, qAlignSize, qstrand, qScaffSize, blocks, mismatch) = line.split()[:13]

                alignedRef[scaff_reference].append((int(refstart),int(refstart)+int(refAlignSize)-1))
                alignedQuery[scaff_query].append((int(qstart),int(qstart)+int(qAlignSize)-1,qstrand))

    return alignedRef,alignedQuery

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


def createtmpBedFiles(alnRef,alntasnformedQuery,referenceGenomeTable, alignedGenomeTable):
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
    iterable_query = b.merge()#.genome_coverage(g=alignedGenomeTable)

    #with open("test3.txt",'w') as out:
    #    for k in iterable_ref:
    #        out.write(str(k))
    return iterable_ref,iterable_query

def generateStats(it_ref,it_query,refDict, queryDict,refGenTable,querGenTable,forwdAlnScaf,revAlnScaf,mixAlnScaf, outBasename):
    with open(outBasename,'w') as outfile:
        for line in it_ref:
            outfile.write(str(line))

def main():
    parser = argparse.ArgumentParser(description='Script to analyse lastal tab output file.')
    parser.add_argument(dest='lastOutput', metavar='lastTabOut', help='Lastal output file in TAB format mapping the reference inputted to the query.')
    parser.add_argument(dest='referenceGenomeTable', metavar='refGenomeTable', help='Tab delimited file displaying reference sequences and their length.')
    parser.add_argument(dest='alignedGenomeTable', metavar='alnGenomeTable',  help='Tab delimited file dsiplaying query sequences and their length.')
    parser.add_argument(dest='outputBasemame', metavar='outputBasename', help='Basename to write output files.')
    args = parser.parse_args()

    dictRef,dictQuery=createGenomeTableDict(args.referenceGenomeTable, args.alignedGenomeTable)
    print("Processing lastal output..")
    alignRef,alignQuery=processLastal(args.lastOutput,dictRef,dictQuery)
    print("Transforming negative strand alignments")
    transformedQueryDic,totallyAlnForwd,totallyAlnRev,mixedAlnFraction= transformNegativeStrandCoordinates(alignQuery,dictQuery)
    print("Creating tmp bed files ..")
    it_ref,it_query = createtmpBedFiles(alignRef,transformedQueryDic,args.referenceGenomeTable, args.alignedGenomeTable)
    print("Generating stats..")
    generateStats(it_ref,it_query,alignRef,transformedQueryDic,dictRef,dictQuery,totallyAlnForwd,totallyAlnRev,mixedAlnFraction,args.outputBasemame)



if __name__ == "__main__":
    main()
