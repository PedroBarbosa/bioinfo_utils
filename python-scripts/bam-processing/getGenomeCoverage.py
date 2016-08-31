__author__ = 'pedro'

import argparse
import subprocess
import logging
import sys
import os
import csv
from collections import defaultdict
from collections import OrderedDict
import numpy
from tempfile import NamedTemporaryFile
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')



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

#            teste=[]
#            for i in bitflagList:
#                if not i in teste:
#                    teste.append(i)

#            for i in teste:
#                print(i)


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

#    teste=[]
#    for i in bitflagList:
#        if not i in teste:
#            teste.append(i)

#    for i in teste:
#        print(i)
#    return bitflagList
    return tmpfile_bitFlag

def replaceBedScoreToBitFlag(bedFile,bitflag_tmp):

    out = os.path.abspath(bedFile) + "BitFlag"
    with open(out,'w+') as out_file:

        sub=subprocess.Popen(["awk",'FNR==NR{a[NR]=$0;next}{$5=a[FNR]};OFS="\t"',bitflag_tmp, bedFile],stderr = subprocess.PIPE, stdout = out_file)
        sub.wait()
        subprocess.Popen(["mv", out, bedFile])


def genomeCoverageFromBed(bedfile,genome,stranded):


    with open(genome) as fin:
        genome_size = sum(int(r[1]) for r in csv.reader(fin, delimiter = "\t"))
    fin.close()

    subTotal = subprocess.check_output(["bedtools", "genomecov", "-i", bedfile, "-g", genome])
    stdout = subTotal.decode("utf-8").split("\n")
    logging.info("\tGenerating stats for coverage..")
    generateStatsAndCoverageFiles(stdout,bedfile,-1,genome_size)

    if stranded:
        logging.info("\tCalculating genome coverage in each strand independently..")
        for i in range(0,1):
            if i == 0:
                subTotal = subprocess.check_output(["bedtools", "genomecov", "-i", bedfile, "-g", genome, "-strand", "+"])
                stdout = subTotal.decode("utf-8").split("\n")
                logging.info("\tGenerating stats for coverage in the positive strand..")
                generateStatsAndCoverageFiles(stdout,bedfile,i,genome_size)
            elif i == 1:
                subTotal = subprocess.check_output(["bedtools", "genomecov", "-i", bedfile, "-g", genome, "-strand", "-"])
                stdout = subTotal.decode("utf-8").split("\n")
                logging.info("\tGenerating stats for coverage in the  negative strand..")
                generateStatsAndCoverageFiles(stdout,bedfile,i,genome_size)



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
                print(attributes)
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
                        print("estou cá, espero que só 1 vez")
                        del dictGenomeCoverageBedtools[previous_scaff][-4:]
                        dictGenomeCoverageBedtools[previous_scaff].extend([round(regionsWithCov,4),round(cov5,4),round(cov10,4),round(cov50,4)])


                    #if last scaffold had 0 of coverage
                    #else:
                    #    print("último scaffold,não tenho cobertura nenhuma")
                    #    dictGenomeCoverageBedtools[previous_scaff].extend([0,0,0,0])

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


    logging.info("Genome size\t%i" % genomeSize)
    logging.info("Genome fraction with coverage 0\t%s" % genomeCov0)
    for k,v in iter(dictWholeGenomeCovLarger0.items()):
        logging.info("Genome fraction with coverage %s\t%s" % (k,v))

    writeDict(dictGenomeCoverageBedtools, bedfile,i)
#    for k,v in iter(dictGenomeCoverageBedtools.items()):
#        if v[1] != '1':
#            print(k,v)
#    print(len(dictGenomeCoverageBedtools))


def writeDict(dict,bedfile, i):
    out_file = ""
    if i == 0:
        out_file = os.path.join(os.path.dirname(os.path.realpath(bedfile)), os.path.basename(bedfile).split('.bed')[0] + '-indScaffCov-plusStrand.txt')
    elif i == 1:
        out_file = os.path.join(os.path.dirname(os.path.realpath(bedfile)), os.path.basename(bedfile).split('.bed')[0] + '-indScaffCov-minusStrand.txt' )
    else:
        out_file = os.path.join(os.path.dirname(os.path.realpath(bedfile)), os.path.basename(bedfile).split('.bed')[0] + '-indScaffCov-bothStrands.txt' )

    with open(out_file, 'w+') as fileout:
        fileout.write('\t'.join(['#scaffold_id','#scaffold_length','#fraction 0 cov', '#fraction > 0 cov', '#fraction 5 > cov < 10', '#fraction 10 > cov < 50', '#fraction > 50 cov']))

        for scaffold, cov in sorted(dict.items(), key=lambda x: x[1][2], reverse=True):
            fileout.write('\n' + scaffold + '\t' + '\t'.join(str(c) for c in cov))


class AlignmentsAnalysis(object):
    def __init__(self, bedFile):
        self.bedfile = bedFile
def main():

    parser = argparse.ArgumentParser(description='Script to analyse the genome coverage and orientation of the scaffolds. Requires samtools and bedtools to be on the system path. BAM files must be sorted.')
    parser.add_argument(dest='bam_files', metavar='bamFiles', nargs='+', help='Bam files to process.')
    parser.add_argument('-m', metavar='mapper', required = True, choices=['bwa_aln','bwa_mem','star','bowtie2'], help='Mapper used in the alignments. Available choices: [bwa_aln,bwa_mem,star,bowtie2].')
    parser.add_argument('-gn', metavar='genomeTable', required = True, help='Tab delimited genome file. Ex:chrom_name    size(bp). Required by bedtools.')
    parser.add_argument('-o', metavar='outputDirectory', required = True, help='Output directory to write the results')
    parser.add_argument('-n', '--nofilterBam', action='store_true', help='If set, the filtering of BAM files will not be perfomed. Default: Process bam files.' )
    parser.add_argument('-s','--singleEnd', action='store_true', help='Single end read mappings. Default:Paired-end')
    parser.add_argument('-str','--stranded', action='store_true', help='Perform also strand specific coverage analysis. Default: Only both strands analysis.')

    args = parser.parse_args()

    for file in args.bam_files:
        if not os.path.isfile(file):
            logging.error("%s file is not valid." % file)
            exit(1)
        if os.path.splitext(file)[1] != ".bam":
            logging.error("Only files with .bam extension are allowed. Exiting!")
            exit(1)

    if not os.path.exists(args.o):
        os.makedirs(args.o)

    for file in args.bam_files:
        logging.info("Processing %s file.." % os.path.basename(file))
        outBed = os.path.join(args.o,os.path.basename(file).replace('.bam', '.bed'))

        if not args.nofilterBam:
            if args.singleEnd:
                bitflag_tmp = filterSingleEndBamAlignments(file,args.o,args.m)
                logging.info("\tAdding bitFlags to the score column in BED file..")
                replaceBedScoreToBitFlag(outBed,bitflag_tmp.name)
                logging.info("Calculating genome coverage..")
                genomeCoverageFromBed(outBed, args.gn, args.stranded)
            else:
                bitflag_tmp = filterBamAlignments(file,args.o,args.m)
                logging.info("\tAdding bitFlags to the score column in BED file..")
                replaceBedScoreToBitFlag(outBed,bitflag_tmp.name)
                logging.info("Calculating genome coverage..")
                genomeCoverageFromBed(outBed, args.gn,args.stranded)

        else:
            logging.info("No filtering of BAM files will be done.")
            bitflag_tmp = bamToBedFromBam(file,args.o,args.m)
            logging.info("Adding bitFlags to the score column in BED file..")
            replaceBedScoreToBitFlag(outBed,bitflag_tmp.name)
            logging.info("Calculating genome coverage..")
            genomeCoverageFromBed(outBed, args.gn, args.stranded)

        logging.info("DONE!!")

if __name__ == "__main__":
    main()

