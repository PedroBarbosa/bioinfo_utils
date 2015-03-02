import argparse
import sys
import re
import numpy as np

__author__ = 'pedro'


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        self.print_help()
        sys.exit(2)



def extractIdentity(blastFile,outputFile):

    print "Parsing blastP file..\n"
    open(outputFile, 'w').close() ## empty output file

    with open(outputFile, "a") as out:
        out.write("GeneID" + "\t" + "Best_hit_identity"  + "\t" + "Alignment_length" + "\t" + "Alignment_mismatch" + "\n")

    noSimilarities = 0
    noSimilaritiesList, identityList = [],[]
    isnewGene = False
    with open(blastFile) as file:
        for line in file:

            line = line.strip()

            if line.startswith('Query='):
                isnewGene = True # To get identity only for the best hit
                query = line.split('Query=',1)[1][1:]  #get query
                continue


            elif "Identities =" in line and isnewGene:

                isnewGene = False
                alignment_features = line.split(',')
                identity = re.findall("\((\d+).*", alignment_features[0])[0] #get identity

                mismatch = 0
                aln_length = 0

                m = re.search('(\d+)\/(\d+)', alignment_features[0])

                if m:
                    mismatch = int(m.group(2)) - int(m.group(1))
                    aln_length = m.group(2)
                else:
                    mismatch_value = ""
                    aln_length = ""


                identityList.append(int(identity))
                with open(outputFile, "a") as out:
                    out.write(query + "\t" + identity + "%"  + "\t" + aln_length + "\t" + str(mismatch) + "\n")


            elif "No hits found" in line and isnewGene:
                noSimilarities +=1
                noSimilaritiesList.append(query)


    for i in noSimilaritiesList:
        with open(outputFile, "a") as out:
            out.write(i + "\t" + "No_hits" + "\n")

    total_genes = len(noSimilaritiesList) + len(identityList)
    hits_found = len(identityList)
    average_ident = round(np.mean(identityList),3)
    max_ident =  max(identityList)
    min_ident = min(identityList)


    with open(outputFile, "a") as out:
        out.write("\n" + "\n")
        out.write("Total number of genes processed: " + str(total_genes))
        out.write("\nGenes with similarities: " + str(hits_found))
        out.write("\nAverage identity: " + str(average_ident) + "%")
        out.write("\nMaximum identity found: " + str(max_ident) + "%")
        out.write("\nMinimum identity found: " + str(min_ident) + "%")

    print "Total number of genes processed: ", total_genes
    print "Genes with similarities: ", hits_found
    print "Average identity of the best hit for all genes: ", average_ident, "%"
    print "Maximum identity found: ", max_ident,"%"
    print "Minimum identity found: ", min_ident,"%"

parser = MyParser()
parser.add_argument('-i', dest='blastFile',metavar='blast file', required=True, help='Original file with the blast results.')
parser.add_argument('-o', dest='outputFile', metavar='output file', required = True, help='Output file to write the results.')
args = parser.parse_args()



if __name__ == "__main__":

    extractIdentity(args.blastFile, args.outputFile)

