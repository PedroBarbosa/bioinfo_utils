__author__ = 'pedro'

import argparse
import subprocess
import logging
import sys
import os
import csv
from collections import defaultdict
from collections import OrderedDict
from operator import itemgetter
from tempfile import NamedTemporaryFile
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')


global dictGenomeCoverageBedtools
def filterBamAlignments(bamfile,outDir,mapper):
    logging.info("Filtering %s file using the specific %s flags and converting to BED format.." % (os.path.basename(bamfile),mapper))
    tmpfile = NamedTemporaryFile(delete=True)
    tmpfile_bitFlag = NamedTemporaryFile(delete=True)
    if mapper == 'bwa_aln':
        logging.info("\tFiltering out secondary alignments and low quality mappings..")
        UMR0 = subprocess.Popen(["samtools","view","-h","-F260",bamfile], stdout=subprocess.PIPE)
        UMR1 = subprocess.Popen(["grep","-w","X0:i:1"], stdin=UMR0.stdout, stdout=subprocess.PIPE)

        with open(tmpfile.name,'w+') as temp_file:
            logging.info("\tFiltering by specific mapper flags and creating tmp file of the result..")
            UMR = subprocess.Popen(["grep","-w","X1:i:0"], stdin=UMR1.stdout, stdout=temp_file)
            UMR.wait()
            UMRbitflag = subprocess.Popen(["samtools", "view", "-S", tmpfile.name], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            UMRbam = subprocess.Popen(["samtools", "view", "-Shb", tmpfile.name], stderr=subprocess.PIPE, stdout=subprocess.PIPE)

            logging.info("\tExtracting bitflags ..")
            alignments = UMRbitflag.stdout.readlines()
            with open(tmpfile_bitFlag.name, 'w+') as tmp_file_bitflag:
                for line in alignments:
                    decode = line.decode("utf-8")
                    tmp_file_bitflag.write(decode.split("\t")[1] + '\n')

        tmp_file_bitflag.close()
        logging.info("\tConverting bam to bed using BEDtools..")
        with open(os.path.join(outDir,os.path.basename(bamfile).replace('.bam', '.bed')), 'w+') as out_file:
            bam2bed = subprocess.Popen(["bedtools", "bamtobed", "-cigar", "-i", "stdin"], stdin=UMRbam.stdout, stdout=out_file)

            bam2bed.wait()
            out_file.seek(0)
            logging.info("DONE!!! %i alignments kept for further analysis." % len(out_file.readlines()))

            out_file.close()
            UMR0.stdout.close()
            UMR1.stdout.close()


    elif mapper == 'bwa_mem':
        logging.info("\tFiltering out secondary alignments and low quality mappings..")
        q10 = subprocess.Popen(["samtools","view", "-h", "-F260", "-q10", bamfile], stdout=subprocess.PIPE)
        with open(tmpfile.name,'w+') as temp_file:
            logging.info("\tFiltering by specific mapper flags and creating tmp file of the result..")
            q10XA = subprocess.Popen(["grep", "-v", "XA:"], stdin=q10.stdout, stdout=temp_file)
            q10XA.wait()

            q10XAbitflag = subprocess.Popen(["samtools", "view", "-S", tmpfile.name], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            q10XAbam = subprocess.Popen(["samtools", "view", "-Shb", tmpfile.name], stderr=subprocess.PIPE, stdout=subprocess.PIPE)

            logging.info("\tExtracting bitflags ..")
            alignments = q10XAbitflag.stdout.readlines()
            with open(tmpfile_bitFlag.name, 'w+') as tmp_file_bitflag:
                for line in alignments:
                    decode = line.decode("utf-8")
                    tmp_file_bitflag.write(decode.split("\t")[1] + '\n')
        tmp_file_bitflag.close()
        logging.info("\tConverting bam to bed using BEDtools..")
        with open(os.path.join(outDir,os.path.basename(bamfile).replace('.bam', '.bed')), 'w+') as out_file:
            bam2bed = subprocess.Popen(["bedtools", "bamtobed", "-cigar", "-i", "stdin"], stdin=q10XAbam.stdout, stdout=out_file)

            bam2bed.wait()
            out_file.seek(0)
            logging.info("DONE!!! %i alignments kept for further analysis." % len(out_file.readlines()))

            out_file.close()
            q10.stdout.close()
            q10XAbam.stdout.close()


    elif mapper == 'bowtie2':
        logging.info("\tFiltering out secondary alignments and low quality mappings..")
        q10 = subprocess.Popen(["samtools", "view", "-h", "-F260", "-q10", bamfile],stdout=subprocess.PIPE)

        with open(tmpfile.name,'w+') as temp_file:
            logging.info("\tFiltering by specific mapper flags and creating tmp file of the result..")
            q10XS = subprocess.Popen(["grep", "-v", "XS:"], stdin=q10.stdout, stdout=temp_file)
            q10XS.wait()

            q10XSbitflag = subprocess.Popen(["samtools", "view", "-S", tmpfile.name], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            logging.info("\tExtracting bitflags..")
            alignments = q10XSbitflag.stdout.readlines()
            with open(tmpfile_bitFlag.name, 'w+') as tmp_file_bitflag:
                for line in alignments:
                    decode = line.decode("utf-8")
                    tmp_file_bitflag.write(decode.split("\t")[1] + '\n')
        tmp_file_bitflag.close()
        logging.info("\tConverting bam to bed using BEDtools..")
        q10XSbam = subprocess.Popen(["samtools", "view", "-Shb", tmpfile.name], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        with open(os.path.join(outDir,os.path.basename(bamfile).replace('.bam', '.bed')), 'w+') as out_file:
            bam2bed = subprocess.Popen(["bedtools", "bamtobed", "-cigar", "-i", "stdin"], stdin=q10XSbam.stdout, stdout=out_file)
            bam2bed.wait()
            out_file.seek(0)
            logging.info("DONE!!! %i alignments used for further analysis." % len(out_file.readlines()))

            out_file.close()
            q10.stdout.close()
            q10XSbam.stdout.close()
            q10XSbitflag.stdout.close()


    elif mapper == 'star':

        with open(tmpfile.name,'w+') as temp_file:
            logging.info("\tExtracting unique mapped reads from STAR alignments to tmp file..")
            q255 = subprocess.Popen(["samtools","view","-h", "-q255", bamfile], stdout=temp_file)
            q255.wait()


            q255bitflag = subprocess.Popen(["samtools", "view", "-S", tmpfile.name], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            q255bam = subprocess.Popen(["samtools", "view", "-Shb", tmpfile.name], stderr=subprocess.PIPE, stdout=subprocess.PIPE)

            logging.info("\tExtracting bitflags ..")
            alignments = q255bitflag.stdout.readlines()
            with open(tmpfile_bitFlag.name, 'w+') as tmp_file_bitflag:
                for line in alignments:
                    decode = line.decode("utf-8")
                    tmp_file_bitflag.write(decode.split("\t")[1] + '\n')

        tmp_file_bitflag.close()
        logging.info("\tConverting bam to bed using BEDtools..")
        with open(os.path.join(outDir,os.path.basename(bamfile).replace('.bam', '.bed')), 'w+') as out_file:
            bam2bed = subprocess.Popen(["bedtools", "bamtobed", "-split","-cigar", "-i", "stdin"], stdin=q255bam.stdout, stdout=out_file)

            bam2bed.wait()
            out_file.seek(0)
            logging.info("DONE!!! %i alignments kept for further analysis." % len(out_file.readlines()))

        q255bam.stdout.close()
        out_file.close()

    return tmpfile_bitFlag

def filterSingleEndBamAlignments(bamfile, outDir, mapper):

    logging.info("Filtering %s file using the specific %s flags and converting to BED format.." % (os.path.basename(bamfile),mapper))
    tmpfile = NamedTemporaryFile(delete=True)
    tmpfile_bitFlag = NamedTemporaryFile(delete=True)

    if mapper == 'bowtie2':

        with open(tmpfile.name,'w+') as temp_file:

            logging.info("\tCreating tmp file of the filtered alignments..")
            q255 = subprocess.Popen(["samtools","view","-h", "-q255", bamfile], stdout=temp_file)

            logging.info("\tExtracting bitflags to list..")
            q255bitflag = subprocess.Popen(["samtools","view","-S", "-q255", tmpfile.name], stdout=subprocess.PIPE)
            alignments = q255bitflag.stdout.readlines()
            with open(tmpfile_bitFlag.name, 'w+') as tmp_file_bitflag:
                for line in alignments:
                    decode = line.decode("utf-8")
                    tmp_file_bitflag.write(decode.split("\t")[1]+ '\n')

        tmp_file_bitflag.close()
        logging.info("\tConverting bam to bed using BEDtools..")
        q255bam = subprocess.Popen(["samtools", "view", "-Shb", tmpfile.name], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        with open(os.path.join(outDir,os.path.basename(bamfile).replace('.bam', '.bed')), 'w+') as out_file:
            bam2bed = subprocess.Popen(["bedtools", "bamtobed","-cigar", "-i", "stdin"], stdin=q255bam.stdout, stdout=out_file)
            bam2bed.wait()
            out_file.seek(0)
            logging.info("DONE!!! %i alignments used for further analysis." % len(out_file.readlines()))


        q255bam.stdout.close()
        out_file.close()

    else:
        logging.error("Processing of the single end reads mapping is only available for bowtie2. Sorry!")
        exit(1)

    return tmpfile_bitFlag


def bamToBedFromBam(bamfile, outputDir, mapper):

    tmpfile_bitFlag = NamedTemporaryFile(delete=True)
    logging.info("Removing just the unmapped reads and extracting bitflags to list.")
    bitflag = subprocess.Popen(["samtools", "view", "-F4", bamfile], stdout=subprocess.PIPE)
    alignments = bitflag.stdout.readlines()

    with open(tmpfile_bitFlag.name, 'w+') as tmp_file_bitflag:
        for line in alignments:
            decode = line.decode("utf-8")
            tmp_file_bitflag.write(decode.split("\t")[1]+'\n')

    tmp_file_bitflag.close()
    logging.info("Converting BAM %s to BED.." % os.path.basename(bamfile))
    with open(os.path.join(outputDir,os.path.basename(bamfile).replace('.bam', '.bed')), 'w+') as out_file:

        if mapper == "star": #if RNA-seq mappings
            bam2bed = subprocess.Popen(["bedtools", "bamtobed", "-cigar", "-split", "-i", bamfile], stdout=out_file)
            bam2bed.wait()
            out_file.seek(0)
            logging.info("DONE!!! %i alignments used for further analysis." % len(out_file.readlines()))
            out_file.close()

        else:
            bam2bed = subprocess.Popen(["bedtools", "bamtobed", "-cigar", "-i", bamfile], stdout=out_file)
            bam2bed.wait()
            out_file.seek(0)
            logging.info("DONE!!! %i alignments used for further analysis." % len(out_file.readlines()))
            out_file.close()

    return tmpfile_bitFlag

def replaceBedScoreToBitFlag(bedFile,bitflag_tmp):

    out = os.path.abspath(bedFile) + "BitFlag"
    with open(out,'w+') as out_file:

        sub=subprocess.Popen(["awk",'FNR==NR{a[NR]=$0;next}{$5=a[FNR]};OFS="\t"',bitflag_tmp, bedFile],stderr = subprocess.PIPE, stdout = out_file)
        sub.wait()
        subprocess.Popen(["mv", out, bedFile])


def genomeCoverageFromBed(bedfile,genome,stranded,outdir,alnType,noCovAnalysis,noAlignmentAnalysis):
    with open(genome) as fin:
        genome_size = sum(int(r[1]) for r in csv.reader(fin, delimiter = "\t"))
    fin.close()
    dicCov={}
    if type(bedfile) is list:#when input is bed
        for file in bedfile:
            if noCovAnalysis == False:
                logging.info("[Genome coverage analysis started]: Calculating genome coverage of %s file directly from a bed input." % file)
                dicCov=bedtoolsCMD(file,genome,stranded,genome_size)
            if noAlignmentAnalysis == False:
                ob=AlignmentsAnalysis(file,outdir,alnType,noCovAnalysis,dicCov)
                ob.processGenomeTable(genome)
                ob.bedAnalysis()
    else:#when input is bam
        if noCovAnalysis == False:
            dicCov=bedtoolsCMD(bedfile,genome,stranded,genome_size)
        if noAlignmentAnalysis == False:
            ob=AlignmentsAnalysis(bedfile,outdir,alnType,noCovAnalysis,dicCov)
            ob.processGenomeTable(genome)
            ob.bedAnalysis()

def bedtoolsCMD(file,genome,stranded,genome_size):

    subTotal = subprocess.check_output(["bedtools", "genomecov", "-i", file, "-g", genome])
    stdout = subTotal.decode("utf-8").split("\n")
    logging.info("\tGenerating stats for coverage..")
    dicCov=generateStatsAndCoverageFiles(stdout,file,-1,genome_size)

    if stranded:
        logging.info("\tCalculating genome coverage in each strand independently..")
        for i in range(0,1):
            if i == 0:
                subTotal = subprocess.check_output(["bedtools", "genomecov", "-i", file, "-g", genome, "-strand", "+"])
                stdout = subTotal.decode("utf-8").split("\n")
                logging.info("\tGenerating stats for coverage in the positive strand..")
                generateStatsAndCoverageFiles(stdout,file,i,genome_size)
            elif i == 1:
                subTotal = subprocess.check_output(["bedtools", "genomecov", "-i", file, "-g", genome, "-strand", "-"])
                stdout = subTotal.decode("utf-8").split("\n")
                logging.info("\tGenerating stats for coverage in the  negative strand..")
                generateStatsAndCoverageFiles(stdout,file,i,genome_size)
    return dicCov

def generateStatsAndCoverageFiles(std_out,bedfile,i,genomeSize):
    dictGenomeCoverageBedtools = defaultdict(list)
    dictWholeGenomeCovLarger0 = OrderedDict()
    dictPerScaffoldCovLarger0 = OrderedDict()
    previous_scaff,tuple,genomeCov0,regionsWithCov, cov5, cov10, cov50, attributes = "",(),0,0.0,0.0,0.0,0.0,[]
    for line in std_out:
        attributes = line.split('\t')
        if len(attributes) > 1:
            #whole genome stats
            if attributes[0] == "genome":
                if int(attributes[1]) == 0:
                    genomeCov0 = attributes[4]
                    #Last scaffold:
                    if dictPerScaffoldCovLarger0:
                        for k,v in iter(dictPerScaffoldCovLarger0.items()):
                            if int(k) < 5:
                                regionsWithCov += float(v)
                            elif int(k) >= 5 and int(k) < 10:
                                cov5 += float(v)
                                regionsWithCov += float(v)
                            elif int(k) >= 10 and int(k) < 50:
                                cov10 += v
                                regionsWithCov += float(v)
                            else:
                                cov50 += float(v)
                                regionsWithCov += float(v)

                        del dictGenomeCoverageBedtools[previous_scaff][-4:]
                        dictGenomeCoverageBedtools[previous_scaff].extend([round(regionsWithCov,4),round(cov5,4),round(cov10,4),round(cov50,4)])

                #whole genome non 0 coverage frequencies
                else:
                    dictWholeGenomeCovLarger0[attributes[1]] = attributes[4]

            #new scaffold, report length and region with 0 coverage
            elif not attributes[0] in dictGenomeCoverageBedtools:
                #update previous scaffold with values from non 0 coverage
                if dictPerScaffoldCovLarger0:

                    for k,v in iter(dictPerScaffoldCovLarger0.items()):
                        if int(k) < 5:
                            regionsWithCov += float(v)
                        elif int(k) >= 5 and int(k) < 10:
                            cov5 += float(v)
                            regionsWithCov += float(v)
                        elif int(k) >= 10 and int(k) < 50:
                            cov10 += float(v)
                            regionsWithCov += float(v)
                        else:
                            cov50 += float(v)
                            regionsWithCov += float(v)

                    del dictGenomeCoverageBedtools[previous_scaff][-4:]
                    dictGenomeCoverageBedtools[previous_scaff].extend([round(regionsWithCov,4),round(cov5,4),round(cov10,4),round(cov50,4)])
                    dictPerScaffoldCovLarger0 = OrderedDict()
                    regionsWithCov = 0
                    cov5 = 0
                    cov10 = 0
                    cov50 = 0
                #if new scaffold, report length and region with 0 coverage
                dictGenomeCoverageBedtools[attributes[0]].extend([attributes[3],float(attributes[4]),0,0,0,0])
            #update dict with values of covarage and proportion
            else:
                previous_scaff = attributes[0]
                dictPerScaffoldCovLarger0[attributes[1]] = attributes[4]

    logging.info("\tFinal: Genome size:\t%i" % genomeSize)
    logging.info("\tFinal: Genome fraction with coverage 0:\t%s" % genomeCov0)
    logging.info("\tFinal: Genome fraction with coverage between 1 and 5:\t%s" % sum([float(v) for k,v in dictWholeGenomeCovLarger0.items() if int(k) >=1 and int(k) <=5 ]))
    logging.info("\tFinal: Genome fraction with coverage between 5 and 10:\t%s" % sum([float(v) for k,v in dictWholeGenomeCovLarger0.items() if int(k) >5 and int(k) <=10 ]))
    logging.info("\tFinal: Genome fraction with coverage between 10 and 50:\t%s" % sum([float(v) for k,v in dictWholeGenomeCovLarger0.items() if int(k) >10 and int(k) <=50 ]))
    logging.info("\tFinal: Genome fraction with coverage larger than 50:\t%s" % sum([float(v) for k,v in dictWholeGenomeCovLarger0.items() if int(k) >50 ]))
    #for k,v in iter(dictWholeGenomeCovLarger0.items()):
    #    logging.info("Genome fraction with coverage %s\t%s" % (k,v))
    writeDict(dictGenomeCoverageBedtools, bedfile,i)
    return dictGenomeCoverageBedtools


def writeDict(dict,bedfile, i):
    out_file = ""
    if i == 0:
        out_file = os.path.join(os.path.dirname(os.path.realpath(bedfile)), os.path.basename(bedfile).split('.bed')[0] + '-indScaffCov-plusStrand.txt')
    elif i == 1:
        out_file = os.path.join(os.path.dirname(os.path.realpath(bedfile)), os.path.basename(bedfile).split('.bed')[0] + '-indScaffCov-minusStrand.txt' )
    else:
        out_file = os.path.join(os.path.dirname(os.path.realpath(bedfile)), os.path.basename(bedfile).split('.bed')[0] + '-indScaffCov-bothStrands.txt' )

    with open(out_file, 'w+') as fileout:
        fileout.write('\t'.join(['#scaffold_id','#scaffold_length','#fraction_0_cov', '#fraction_>_0 cov', '#fraction 5_>_cov_<_10', '#fraction_10_>_cov_<_50', '#fraction_>_50_cov']))

        for scaffold, cov in sorted(dict.items(), key=lambda x: x[1][2], reverse=True):
            fileout.write('\n' + scaffold + '\t' + '\t'.join(str(c) for c in cov))


class AlignmentsAnalysis(object):
    def __init__(self, bedFiles,outdir,alnType,noCovAnalysis,dicCov):
        self.bedfiles = bedFiles
        self.noCovAnalysis=noCovAnalysis
        self.gnTableDict = {}
        self.safeBitFlagsFR = {99:'Proper FR alignment 1st pair', 147:'Proper FR alignment 2nd pair',}
        self.safeBitFlagsRF = {83:'Proper RF alignment 1st pair', 163:'Proper RF alignment 2nd pair'}
        self.wrongInsertDistProperOrientation = {161:'Paired RF alignment 2nd pair',81:'Paired RF alignment 1st pair',97:'Paired FR alignment 1st pair',
                                           145:'Paired FR alignment 2nd pair'}
        self.wrongInsertDistWrongOrientation = {65:'Paired FF alignment 1st pair', 129:'Paired FF alignment 2nd pair', 113:'Paired RR alignment 1st pair',
                                                177:'Paired RR alignment 1st pair'}
        self.unmappedPair = {73:'1st pair forward, mate unmapped', 137:'2nd pair forward, mate unmapped', 89:'1st pair reverse,mate unmapped',
                             153:'2nd pair reverse, mate unmapped'}
        self.unmappedPairStrangeBehavior = {121:'1st pair reverse, mate unmapped reverse', 185:'2nd pair reverse, mate unmapped reverse'}
        self.chimeric = {2193:'Paired FR, 2nd pair, supplementary alignment', 2145:'Paired FR, 1st pair, supplementary alignment', 2129:'Paired RF, 1st pair, supplementary alignment',
                         2163:'Proper RR, 1st pair, supplementary alignment', 2113:'Paired FF, 1st pair, supplementary alignment', 2177:'Paired FF, 2nd pair, supplementary alignment',
                         2225:'Paired RR, 2nd pair, supplementary alignment', 2209:'Paired RF, 2nd pair, supplementary alignment',2115:'Proper FF,1st pair, supplementary alignment',
                         2179:'Proper FF, 2nd pair, supplementary alignment', 2195:'Proper FR, 2nd pair,supplementary alignment',2147:'Proper FR, 1st pair, supplementary alignment',
                         2131:'Proper RF, 1st pair, supplementart alignment',2227:'Proper RR, 2nd pair, supplementart alignment', 2211:'RF, 2nd pair, supplementary alignment',
                         2161:'Paired RR, 1st pair, supplementaru alignment',2233:'Mate umapped,RR pair,',2185:'Mate unmapped, 2nd pair forward, supplementary alignment',
                         }
        self.chimericStrangeBehavior = {2233:'2nd pair reverse, mate unmapped reverse, supplementary alignment',2169:'1st pair reverse, mate unmapped reverse, supplementary alignment '}
        self.secondary = {}
        self.outDir=outdir
        self.alnType = alnType
        self.finalDictCov=dicCov
        if self.alnType == "pe":
            self.primary=self.safeBitFlagsFR
        else:
            self.primary=self.safeBitFlagsRF
        self.totalAln, self.totalCorrect, self.totalCorrect_Orientation,self.totalChim,self.totalSecondary, self.totalClipped,self.totalClippedWithin, \
        self.totalClippedScaffoldEnds, self.totalClippedNotChimeric,self.totalClippedChimeric = 0,0,0,0,0,0,0,0,0,0

    def processGenomeTable(self,gnFile):
        with open(gnFile, 'r') as infile:
            for line in infile:
                if not line.startswith("#") and not len(line.split("\t")) == 2:
                    logging.error("Please verify genome table file. It should refer to a 2 column tab separated file with contig/scaffold id and its length")
                    exit(1)
                else:
                    id = line.split("\t")[0]
                    length = line.split("\t")[1]
                    self.gnTableDict[id] = int(length)


    def bedAnalysis(self):
        if type(self.bedfiles) is list:
            for bed in self.bedfiles:
                self.deepAlignmentAnalaysis(bed)
        else:
            self.deepAlignmentAnalaysis(self.bedfiles)

    def deepAlignmentAnalaysis(self,bed):
        hash,conflictAlignReadInfo=defaultdict(list),defaultdict(list)
        logging.info("[Alignments analysis started]: Hashing %s bed file.." % bed)
        with open(bed,'r') as infile:
            for line in infile:
                attr=line.rstrip().split()
                if int(attr[4]) in self.wrongInsertDistProperOrientation or int(attr[4]) in self.wrongInsertDistWrongOrientation or int(attr[4]) in self.chimeric:
                    unique_readID=attr[3].split("/")
                    conflictAlignReadInfo[unique_readID[0]].append((attr[0],unique_readID[1],attr[4],attr[6]))

                hash[attr[0]].append((int(attr[1]),int(attr[2]),attr[3],int(attr[4]), attr[6]))
        infile.close()
        logging.info("Processing data..")
        outfile=os.path.join(self.outDir, os.path.basename(bed).replace('.bed','_alignmentInfo.txt'))
        outchim=os.path.join(self.outDir,os.path.basename(bed).replace('.bed','_chimericInfo.txt'))
        outclipping=os.path.join(self.outDir,os.path.basename(bed).replace('.bed','_clippingInfo.txt'))
        with open(outfile, 'w') as out:
            with open(outchim,'w') as outchim:
                with open(outclipping,'w') as outclip:
                    outchim.write("#readID\tmappedScaffold\tpair\tbitflag\tCIGAR\n")
                    outclip.write("#Clipped_Alignments - alignment that displays soft/hard clipping.\n")
                    outclip.write("#Clipped_Chimeric - represents a supplementary alignment of a chimeric one.\n")
                    outclip.write("#Clipped_NonChimeric - represents an alignment that is either the representative record of a chimeric alignment or just a "
                                  "linear alignment that was soft/hard clipped.\n")
                    outclip.write("#Clipped_ScaffEnd - alignment [chimeric or not] that is clipped in the scaffold boundaries.\n")
                    outclip.write("#Clipped_WithinScaff - alignment [chimeric or not] that was clipped within scaffold boundaries.\n")
                    if self.noCovAnalysis:
                        outclip.write("#scaffoldID\tLength\tNumber_alignments\tClipped_Alignments\tClipped_Chimeric\tClipped_NonChimeric\tClipped_ScaffEnd\tClippedWithinScaff\n")
                    else:
                        outclip.write("#scaffoldID\tLength\tFractionCoverage>0\tNumber_alignments\tClipped_Alignments\tClipped_Chimeric\tClipped_NonChimeric\tClipped_ScaffEnd\tClippedWithinScaff\n")
                    out.write("#Number_alignments - total number of alignments found in bam/bed files.\n")
                    out.write("#Correct - pairs align concordantly in the expected insert distance.\n")
                    out.write("#Correct&Orientation - alignments belonging to #Correct, that obbey to the expected corret orientation "
                              "provided by the alignment type: mp or pe.\n")
                    out.write("#Potential_correct - pairs map in different scaffolds with good orientation [FR or RF] but unknown distance.\n" )
                    out.write("#Potential_correct&unexpectedOrientation - pairs map in different scaffolds with wrong orientation that might"
                              "be justified by the wrong reference orientation.\n")
                    out.write("#Potential_correct&unmappedMate - alignments that map perfectly [regardless proper orientation given alignment type],"
                              " but the mate didn't map.\n")
                    out.write("#Wrong_distance - pairs map in the same scaffold with wrong insert size distribution.\n")
                    out.write("#Wrong_orientation - pairs map in the same scaffold with wrong orientation [FF,RR].\n")
                    out.write("#Unknown_fate - pairs for which one mate probably didn't pass mapping quality filters.\n")
                    out.write("#Chimeric - pair fow which a single linear alignment could not be found.\n")
                    if self.noCovAnalysis:
                        out.write("#scaffoldID\tLength\tNumber_alignments\tCorrect [fraction%]\tCorrect&Orientation\tPotential_correct\t"
                                  "Potential_correct&unexpectedOrientation\tPotential_correct&unmappedMate\tWrong_distance\tWrong_orientation\tUnknown_fate\tChimeric\n")
                    else:
                        out.write("#scaffoldID\tLength\tFractionCoverage>0\tNumber_alignments\tCorrect [fraction%]\tCorrect&Orientation\tPotential_correct\t"
                                  "Potential_correct&unexpectedOrientation\tPotential_correct&unmappedMate\tWrong_distance\tWrong_orientation\tUnknown_fate\tChimeric\n")

                    for k,alignments in hash.items():
                        numb_align=len(alignments)
                        self.totalAln+=numb_align
                        corrected, corrected_orientation,potential_correct,potential_correct_wrongOrientation,potential_correct_unmappedMate,wrong_distance, \
                        wrong_orientation,unknown_fate,chimeric,clipped,clippedScafEnd,clippedWithin,clippedChim,clippedNonChim=0,0,0,0,0,0,0,0,0,0,0,0,0,0
                        for i in alignments:
                            if "H" in i[4] or "S" in i[4]:
                                self.totalClipped+=1
                                clipped+=1
                                if i[0] != 0 and i[1] != self.gnTableDict[k]:
                                    self.totalClippedWithin+=1
                                    clippedWithin+=1
                                else:
                                    self.totalClippedScaffoldEnds+=1
                                    clippedScafEnd+=1
                                if not i[3] in self.chimeric and not i[3] in self.chimericStrangeBehavior:
                                    self.totalClippedNotChimeric+=1
                                    clippedNonChim+=1
                                else:
                                    self.totalClippedChimeric+=1
                                    clippedChim+=1

                            if i[3] in self.safeBitFlagsFR or i[3] in self.safeBitFlagsRF:#perfect alignments
                                self.totalCorrect+=1
                                corrected+=1
                                if i[3] in self.primary.keys():
                                    self.totalCorrect_Orientation+=1
                                    corrected_orientation+=1
                            elif i[3] in self.wrongInsertDistProperOrientation or i[3] in self.wrongInsertDistWrongOrientation:#problematic
                                if len(conflictAlignReadInfo[i[2].split("/")[0]]) == 1:#if both pairs aligned,but one of the pairs didn't pass quality filters and was removed, although the bitflag shows the opposite
                                    unknown_fate+=0
                                elif len(conflictAlignReadInfo[i[2].split("/")[0]]) == 2:#if both pairs are uniquely represented
                                    s = set(e[0] for e in conflictAlignReadInfo[i[2].split("/")[0]])

                                    if len(s) == 1 and i[3] in self.wrongInsertDistProperOrientation:#if pairs map in the same scaffold but somehow have wrong distance
                                        wrong_distance+=1
                                    elif len(s) == 1 and i[3] in self.wrongInsertDistWrongOrientation:#if pairs map in the same scaffold but with wrong orientation
                                        wrong_orientation+=1
                                    elif i[3] in self.wrongInsertDistProperOrientation:#pairs map in different scaffolds with the correct expected orientation
                                        potential_correct+=1
                                    else:#pairs map in different scaffolds with a wrong orientation
                                        potential_correct_wrongOrientation+=1

                                else:
                                    if i[2].split("/")[0] in conflictAlignReadInfo:
                                        logging.warning("%s read has at least one secondary/chimeric alignment." % i[2].split("/")[0])

                            elif i[3] in self.unmappedPair or i[3] in self.unmappedPairStrangeBehavior:#unmapped mate
                                potential_correct_unmappedMate+=1

                            elif i[3] in self.chimeric or i[3] in self.chimericStrangeBehavior:
                                self.totalChim+=1
                                chimeric+=1
                                chiminfo=conflictAlignReadInfo[i[2].split("/")[0]]

                                for record in sorted(chiminfo, key=itemgetter(1)):
                                    outchim.write(i[2].split("/")[0]  + '\t' + '\t'.join(record) + '\n')

                            else:
                                print("New FLAG %i. talk with pedro to fix" % i[3])


                        if self.noCovAnalysis:
                            out.write("%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n" % (k,self.gnTableDict[k],numb_align,corrected,corrected_orientation
                                                                                            ,potential_correct,potential_correct_wrongOrientation,
                                                                                            potential_correct_unmappedMate,wrong_distance,wrong_orientation,unknown_fate,chimeric))
                            outclip.write("%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n" % (k,self.gnTableDict[k],numb_align,clipped,clippedChim,clippedNonChim,clippedScafEnd,clippedWithin))
                        else:
                            out.write("%s\t%i\t%.3f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n" % (k,self.gnTableDict[k],self.finalDictCov[k][2],numb_align,corrected,corrected_orientation
                                                                                                ,potential_correct,potential_correct_wrongOrientation,
                                                                                                potential_correct_unmappedMate,wrong_distance,wrong_orientation,unknown_fate,chimeric))
                            outclip.write("%s\t%i\t%.3f\t%i\t%i\t%i\t%i\t%i\t%i\n" % (k,self.gnTableDict[k],self.finalDictCov[k][2],numb_align,clipped,clippedChim,clippedNonChim,clippedScafEnd,clippedWithin))

                    for scaff,length in self.gnTableDict.items():
                        if scaff not in hash:
                            out.write("%s\t%i\t%i\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (scaff,length,0,'*','*','*','*','*','*','*','*','*'))
                            outclip.write("%s\t%i\t%i\t%s\t%s\t%s\t%s\t%s\n" % (scaff,length,0,'*','*','*','*','*'))

                    logging.info("\tFinal: Total number of alignments:\t%i" % self.totalAln)
                    logging.info("\tFinal: Number of correct alignments:\t%i [%.2f%%]" % (self.totalCorrect,round(self.totalCorrect/self.totalAln*100,2)))
                    logging.info("\tFinal: Number of correct alignments given the expected library orientation:\t%i [%.2f%%]" % (self.totalCorrect_Orientation,round(self.totalCorrect_Orientation/self.totalAln*100,2)))
                    logging.info("\tFinal: Number of chimeric alignments [supplementary alignments, representative not displayed in the bitflag field]:"
                                 "\t%i" % self.totalChim)
                    logging.info("\tFinal: Number of secondary alignments:\t%i" % self.totalSecondary)
                    logging.info("\tFinal: Number of alignments with clipped bases:\t%i" % self.totalClipped)
                    logging.info("\tFinal: Number of chimeric alignments [supplementary alignments - 2048 flag] with clipped bases [redundant field]:\t%i" % self.totalClippedChimeric)
                    logging.info("\tFinal: Number of alignments with clipped bases that is not chimeric [20148 flag]:\t%i" % self.totalClippedNotChimeric)
                    logging.info("\tFinal: Number of alignments with clipped bases in the scaffold ends:\t%i" % self.totalClippedScaffoldEnds)
                    logging.info("\tFinal: Number of aligments with clipped bases within scaffold boundaries:\t%i\n\n\n" % self.totalClippedWithin)

            out.close()
            outchim.close()
            outclip.close()

def main():

    parser = argparse.ArgumentParser(description='Script to analyse the genome coverage and orientation of the scaffolds. Requires samtools and bedtools to be on the system path. BAM files must be sorted.')
    parser.add_argument(dest='bam_files', metavar='bamFiles', nargs='+', help='Bam files to process.')
    parser.add_argument('-m', metavar='mapper', required = True, choices=['bwa_aln','bwa_mem','star','bowtie2'], help='Mapper used in the alignments. Available choices: [bwa_aln,bwa_mem,star,bowtie2].')
    parser.add_argument('-gn', metavar='genomeTable', required = True, help='Tab delimited genome file. Ex:chrom_name    size(bp). Required by bedtools.')
    parser.add_argument('-o', metavar='outputDirectory', required = True, help='Output directory to write the results')
    parser.add_argument('-n', '--nofilterBam', action='store_true', help='If set, the filtering of BAM files will not be perfomed. Default: Process bam files.' )
    parser.add_argument('-s','--singleEnd', action='store_true', help='Single end read mappings. Default:Paired-end')
    parser.add_argument('-str','--stranded', action='store_true', help='Perform also strand specific coverage analysis. Default: Only both strands analysis.')
    parser.add_argument('-t','--alnType', required = True, choices=['pe','mp'], help='Type of alignments. If paired-end expects FR orientation, else RF.')
    parser.add_argument('-bed', '--bedFiles', action='store_true', help='Input files are in bed format with the following specification:[ref_id start end readID mapBitflag strand CIGAR].')
    typeOfAnalysis = parser.add_mutually_exclusive_group(required=False)
    typeOfAnalysis.add_argument('-c','--noCoverage', action='store_true', help='Do not execute genome coverage analysis. Stick to the alignments analysis.')
    typeOfAnalysis.add_argument('-a','--noAlignment',action='store_true', help='Do not execute alignments analysis. Stick to the coverage analysis.')
    args = parser.parse_args()

    for file in args.bam_files:
        if not os.path.isfile(file):
            logging.error("%s file is not valid." % file)
            exit(1)
        if not args.bedFiles and os.path.splitext(file)[1] != ".bam":
            logging.error("Only files with .bam extension are allowed. Exiting!")
            exit(1)
        elif args.bedFiles and os.path.splitext(file)[1] != ".bed":
            logging.error("Only files with .bed extension are allowed when '-bed' option is set. Exiting!")
            exit(1)

    if not os.path.exists(args.o):
        os.makedirs(args.o)

    if not args.bedFiles:
        for file in args.bam_files:
            logging.info("Processing %s file.." % os.path.basename(file))
            outBed = os.path.join(args.o,os.path.basename(file).replace('.bam', '.bed'))

            if not args.nofilterBam:
                if args.singleEnd:
                    bitflag_tmp = filterSingleEndBamAlignments(file,args.o,args.m)
                    logging.info("\tAdding bitFlags to the score column in BED file..")
                    replaceBedScoreToBitFlag(outBed,bitflag_tmp.name)
                    logging.info("Calculating genome coverage..")
                    genomeCoverageFromBed(outBed, args.gn, args.stranded,args.o, args.alnType,args.noCoverage,args.noAlignment)
                else:
                    bitflag_tmp = filterBamAlignments(file,args.o,args.m)
                    logging.info("\tAdding bitFlags to the score column in BED file..")
                    replaceBedScoreToBitFlag(outBed,bitflag_tmp.name)
                    logging.info("Calculating genome coverage..")
                    genomeCoverageFromBed(outBed, args.gn,args.stranded,args.o, args.alnType,args.noCoverage,args.noAlignment)

            else:
                logging.info("No filtering of BAM files will be done.")
                bitflag_tmp = bamToBedFromBam(file,args.o,args.m)
                logging.info("Adding bitFlags to the score column in BED file..")
                replaceBedScoreToBitFlag(outBed,bitflag_tmp.name)
                logging.info("Calculating genome coverage..")
                genomeCoverageFromBed(outBed, args.gn, args.stranded,args.o, args.alnType,args.noCoverage,args.noAlignment)

            logging.info("DONE!!")
    else:
        genomeCoverageFromBed(args.bam_files,args.gn,args.stranded,args.o, args.alnType,args.noCoverage,args.noAlignment)


if __name__ == "__main__":
    main()

