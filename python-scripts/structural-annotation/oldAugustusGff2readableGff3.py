import argparse

def createGenomeTableDict(refGenomeTable):
    gnmTbl={}
    with open(refGenomeTable, 'r') as infile:
        for line in infile:
            if len(line.split("\t")) != 2:
                print("Error, genome table files must have exactly 2 columns.")
                exit(1)
            elif line.split("\t")[0] in gnmTbl:
                print("Error, repeated scaffold id in genome table file %s" % refGenomeTable)
            else:
                gnmTbl[line.split("\t")[0]] = line.split("\t")[1].rstrip()
        #print("%i reference sequences processed" % len(gnmTbl))
    infile.close()
    return gnmTbl

def adjustGffFromBedCoordinates(gff,bed,outgff,dictin):
    with open(bed,'r') as bedin:
        dict_bed={}
        for line in bedin:
            if not line.startswith("#"):
                fields = line.rstrip().split("\t")
                if fields[0] in dict_bed:
                    print("#Error: %s sequence ID is present more than once in bed file. Only one trimmed region per sequence is allowed")
                    exit(1)
                else:
                    dict_bed[fields[0]] = [int(fields[1]),int(fields[2])]
    bedin.close()

    with open(gff,'r') as gffin:
        with open(outgff,'w') as outfl:
            outfl.write("##gff-version AUGUSTUS output\n##Output generated with Pedro script to adjust gff features coordinates based on the trimming done on the fasta file in the terminal"
                        " regions of some sequences.\n")
            interval=(0,0)
            toAdjust=False
            geneID=""
            for line in gffin:
                if not line.startswith("#"):
                    fields=line.rstrip().split("\t")
                    if fields[0] not in dict_bed.keys():
                        outfl.write(line)
                    else:
                        if fields[2] == "gene":
                            toAdjust=False
                            geneID=fields[8]
                            gene_interval=(int(fields[3]),int(fields[4]))
                            if dict_bed[fields[0]][0] == 0:
                                if gene_interval[0] < dict_bed[fields[0]][1]: #if beginning of gene is lower than end bed coordinate, we have an issue
                                    print("Serious WARNING: %s gene is located upstream of the region trimmed in the scaffold %s. Perhaps you want to remove "
                                          "this gene??" % (fields[8],fields[0]))
                                else:
                                    toAdjust=True
                                    fields[3] = str(int(fields[3]) - dict_bed[fields[0]][1])
                                    fields[4] = str(int(fields[4]) - dict_bed[fields[0]][1])
                                    outfl.write("\t".join(fields) + "\n")
                                    print("%s gene ID coordinates updated (as well as its child)" % fields[8])
                            elif dict_bed[fields[0]][1] == int(dictin[fields[0]]):
                                if gene_interval[1] > dict_bed[fields[0]][0]: #if end of gene is larger than begiining bed coordinat, we have an issue
                                    print("Serious WARNING: %s gene is located downstream of the region trimmed in the scaffold %s. Perhaps you want to remove "
                                          "this gene??" % (fields[8],fields[0]))


                        elif toAdjust == True and geneID in fields[8]:
                            fields[3] = str(int(fields[3]) - dict_bed[fields[0]][1])
                            fields[4] = str(int(fields[4]) - dict_bed[fields[0]][1])
                            outfl.write("\t".join(fields) + "\n")

def readAugustusGff(infile,outfile,dictin):
    with open(infile,'r') as infl:
        with open(outfile,'w') as outfl:
            outfl.write("##gff-version 3\n##Output generated with Pedro script to convert IDs nomenclature in the 9th attributes field.\n")
            geneID=""
            transcriptID=""
            previous_refid=""
            i=0
            for line in infl:
                if not line.startswith("#"):
                    fields=line.split("\t")
                    refid=fields[0]
                    if refid != previous_refid:
                        if refid in dictin:
                            outfl.write("##sequence-region " + refid + " 1 " + dictin[refid] + "\n")
                        else:
                            print("%s seq id not in genome table. Sequence region will not be present for this reference sequence." % fields[0])

                    if fields[2] == "gene":
                        geneID=fields[8].rstrip()
                        outfl.write('\t'.join(fields[:-1]) + '\tID=' + geneID + "\n")
                    elif fields[2] == "transcript":
                        transcriptID=fields[8].rstrip()
                        fields[2] = "mRNA"
                        outfl.write('\t'.join(fields[:-1]) + '\tID=' + transcriptID + ";Parent=" + geneID + "\n")
                        i=1
                    elif fields[2] == "start_codon":
                        outfl.write('\t'.join(fields[:-1]) + '\tID=' + transcriptID + ":start" + ';Parent=' + transcriptID + "\n")
                    elif fields[2] == "stop_codon":
                        outfl.write('\t'.join(fields[:-1]) + '\tID=' + transcriptID + ":stop" + ';Parent=' + transcriptID + "\n")
                    elif fields[2] == "exon" or fields[2] == "CDS":
                        fields[2] = "exon"
                        outfl.write('\t'.join(fields[:-1]) + '\tID=' + transcriptID + ":exon" + str(i) + ';Parent=' + transcriptID + "\n")
                        i+=1
                    elif fields[2] == "intron":
                        continue
                    previous_refid=refid
    infl.close()
    outfl.close()



def main():
    parser = argparse.ArgumentParser(description='Script to convert native non gff3 augustus format into a readable gff3. If -bed option is set, script'
                                                 'will just perform the adjustment of gff coordinates. In such cases, you need to run again the script'
                                                 'on the update gff to get a new one with the gff3 format')
    parser.add_argument(dest='inputgff', metavar='augustusGff', help='Input file.')
    parser.add_argument(dest='outputfgff', metavar='outgff3', help='Name of the output file.')
    parser.add_argument("-g", metavar='-genomeTable', required=True,help="Genome table file to ouptut ##sequence-region in the gff. Advised by genome validation tools.")
    parser.add_argument("-b", metavar='--bed', help="Bed file representing genome regions that were trimmed from the assembly , whick raises the need to adjust also the gff.")
    args = parser.parse_args()

    outdict = createGenomeTableDict(args.g)
    if args.b:
        adjustGffFromBedCoordinates(args.inputgff,args.b,args.outputfgff,outdict)
    else:
        readAugustusGff(args.inputgff,args.outputfgff,outdict)


if __name__ == "__main__":
    main()
