__author__ = 'pedro'

import argparse
import subprocess
import logging
import sys
import os
import fileinput
from tempfile import NamedTemporaryFile
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')



def filterBamAlignments(bamfile,outDir,mapper):
    logging.info("Filtering %s file using the specific %s flags and converting to BED format.." % (os.path.basename(bamfile),mapper))
    bitflagList = []
    tmpfile = NamedTemporaryFile(delete=True)
    if mapper == 'bwa_aln':
        UMR0 = subprocess.Popen(["samtools","view","-h","-F260",bamfile], stdout=subprocess.PIPE)
        UMR1 = subprocess.Popen(["grep","-w","X0:i:1"], stdin=UMR0.stdout, stdout=subprocess.PIPE)

        with open(tmpfile.name,'w+') as temp_file:
            UMR = subprocess.Popen(["grep","-w","X1:i:0"], stdin=UMR1.stdout, stdout=temp_file)
            UMR.wait()
            UMRbitflag = subprocess.Popen(["samtools", "view", "-S", tmpfile.name], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            UMRbam = subprocess.Popen(["samtools", "view", "-Shb", tmpfile.name], stderr=subprocess.PIPE, stdout=subprocess.PIPE)

            alignments = UMRbitflag.stdout.readlines()
            for line in alignments:
                decode = line.decode("utf-8")
                bitflagList.append(decode.split("\t")[1])



        with open(os.path.join(outDir,os.path.basename(bamfile).replace('.bam', '.bed')), 'w+') as out_file:
            bam2bed = subprocess.Popen(["bedtools", "bamtobed", "-cigar", "-i", "stdin"], stdin=UMRbam.stdout, stdout=out_file)

            bam2bed.wait()
            out_file.seek(0)
            logging.info("DONE!!! %i alignments kept for further analysis." % len(out_file.readlines()))

            out_file.close()
            UMR0.stdout.close()
            UMR1.stdout.close()
            UMR.stdout.close()

    elif mapper == 'bwa_mem':
        q10 = subprocess.Popen(["samtools","view", "-h", "-F260", "-q10", bamfile], stdout=subprocess.PIPE)
        q10XA = subprocess.Popen(["grep", "-v", "XA:"], stdin=q10.stdout, stdout=subprocess.PIPE)
        with open(tmpfile.name,'w+') as temp_file:
            q10XAbam = subprocess.Popen(["samtools", "view", "-Shb", "-"], stdin=q10XA.stdout,  stdout=temp_file)
            q10XAbam.wait()
            q10XAbitflag = subprocess.Popen(["samtools", "view", "-S", tmpfile.name], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            q10XAbam = subprocess.Popen(["samtools", "view", "-Shb", tmpfile.name], stderr=subprocess.PIPE, stdout=subprocess.PIPE)

            alignments = q10XAbitflag.stdout.readlines()
            for line in alignments:
                decode = line.decode("utf-8")
                bitflagList.append(decode.split("\t")[1])

        with open(os.path.join(outDir,os.path.basename(bamfile).replace('.bam', '.bed')), 'w+') as out_file:
            bam2bed = subprocess.Popen(["bedtools", "bamtobed", "-cigar", "-i", "stdin"], stdin=q10XAbam.stdout, stdout=out_file)

            bam2bed.wait()
            out_file.seek(0)
            logging.info("DONE!!! %i alignments kept for further analysis." % len(out_file.readlines()))

            out_file.close()
            q10.stdout.close()
            q10XA.stdout.close()
            q10XAbam.stdout.close()


    elif mapper == 'bowtie2':

        q10 = subprocess.Popen(["samtools", "view", "-h", "-F260", "-q10", bamfile],stdout=subprocess.PIPE)
        with open(tmpfile.name,'w+') as temp_file:

            q10XS = subprocess.Popen(["grep", "-v", "XS:"], stdin=q10.stdout, stdout=temp_file)
            q10XS.wait()

            q10XSbitflag = subprocess.Popen(["samtools", "view", "-S", tmpfile.name], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            q10XSbam = subprocess.Popen(["samtools", "view", "-Shb", tmpfile.name], stderr=subprocess.PIPE, stdout=subprocess.PIPE)


            alignments = q10XSbitflag.stdout.readlines()
            for line in alignments:
                decode = line.decode("utf-8")
                bitflagList.append(decode.split("\t")[1])



        with open(os.path.join(outDir,os.path.basename(bamfile).replace('.bam', '.bed')), 'w+') as out_file:
            bam2bed = subprocess.Popen(["bedtools", "bamtobed", "-cigar", "-i", "stdin"], stdin=q10XSbam.stdout, stdout=out_file)
            bam2bed.wait()
            out_file.seek(0)
            logging.info("DONE!!! %i alignments used for further analysis." % len(out_file.readlines()))

            out_file.close()
            q10.stdout.close()
            q10XSbam.stdout.close()
            q10XSbitflag.stdout.close()

            teste=[]
            for i in bitflagList:
                if not i in teste:
                    teste.append(i)

            for i in teste:
                print(i)
    elif mapper == 'star':



        with open(tmpfile.name,'w+') as temp_file:

            q255 = subprocess.Popen(["samtools","view","-h", "-q255", bamfile], stdout=temp_file)
            q255.wait()

            q255bitflag = subprocess.Popen(["samtools", "view", "-S", tmpfile.name], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            q255bam = subprocess.Popen(["samtools", "view", "-Shb", tmpfile.name], stderr=subprocess.PIPE, stdout=subprocess.PIPE)


            alignments = q255bitflag.stdout.readlines()
            for line in alignments:
                decode = line.decode("utf-8")
                bitflagList.append(decode.split("\t")[1])

        with open(os.path.join(outDir,os.path.basename(bamfile).replace('.bam', '.bed')), 'w+') as out_file:
            bam2bed = subprocess.Popen(["bedtools", "bamtobed", "-split","-cigar", "-i", "stdin"], stdin=q255bam.stdout, stdout=out_file)

            bam2bed.wait()
            out_file.seek(0)
            logging.info("DONE!!! %i alignments kept for further analysis." % len(out_file.readlines()))

        q255.stdout.close()
        q255bam.stdout.close()
        out_file.close()

    return bitflagList

def filterSingleEndBamAlignments(bamfile, outDir, mapper):

    logging.info("Filtering %s file using the specific %s flags and converting to BED format.." % (os.path.basename(bamfile),mapper))
    bitflagList = []
    tmpfile = NamedTemporaryFile(delete=True)

    if mapper == 'bowtie2':

        with open(tmpfile.name,'w+') as temp_file:

            q255 = subprocess.Popen(["samtools","view","-h", "-q255", bamfile], stdout=temp_file)

            q255bitflag = subprocess.Popen(["samtools","view","-S", "-q255", tmpfile.name], stdout=subprocess.PIPE)
            alignments = q255bitflag.stdout.readlines()
            for line in alignments:
                decode = line.decode("utf-8")
                bitflagList.append(decode.split("\t")[1])

        q255bam = subprocess.Popen(["samtools", "view", "-Shb", tmpfile.name], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        with open(os.path.join(outDir,os.path.basename(bamfile).replace('.bam', '.bed')), 'w+') as out_file:
            bam2bed = subprocess.Popen(["bedtools", "bamtobed","-cigar", "-i", "stdin"], stdin=q255bam.stdout, stdout=out_file)
            bam2bed.wait()
            out_file.seek(0)
            logging.info("DONE!!! %i alignments used for further analysis." % len(out_file.readlines()))

        q255.stdout.close()
        q255bam.stdout.close()
        out_file.close()

    else:
        logging.error("Processing of the single end reads mapping is only available for bowtie2. Sorry!")
        exit(1)


    return bitflagList

def bamToBedFromBam(bamfile, outputDir, mapper):
    logging.info("Converting BAM %s to BED.." % os.path.basename(bamfile))
    bitflagList = []

    bitflag = subprocess.Popen(["samtools", "view", "-F4", bamfile], stdout=subprocess.PIPE)
    alignments = bitflag.stdout.readlines()
    for line in alignments:
        decode = line.decode("utf-8")
        bitflagList.append(decode.split("\t")[1])

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

    teste=[]
    for i in bitflagList:
        if not i in teste:
            teste.append(i)

    for i in teste:
        print(i)
    return bitflagList


def replaceBedScoreToBitFlag(bedFile,bitflagList):

    with fileinput.FileInput(bedFile, inplace=True) as bedfile:
        for bitflag,line in zip(bitflagList,bedfile):
            elements = line.split()
            print(line.replace(elements[4], bitflag), end='')

# and 147 - FR
#82 and 163 - RF

def main():

    parser = argparse.ArgumentParser(description='Script to analyse the genome coverage and orientation of the scaffolds. Requires samtools and bedtools to be on the system path. BAM files must be sorted.')
    parser.add_argument(dest='bam_files', metavar='bamFiles', nargs='+', help='Bam files to process.')
    parser.add_argument('-m',metavar = '--mapper', required = True, choices=['bwa_aln','bwa_mem','star','bowtie2'], help='Mapper used in the alignments. Available choices: [bwa_aln,bwa_mem,star,bowtie2].')
    parser.add_argument('-g', metavar = '--genomeTable', required = True, help='Tab delimited genome file. Ex:chrom_name    size(bp). Required by bedtools.')
    parser.add_argument('-o',metavar = '--outputDirectory', required = True, help='Output directory to write the results')
    parser.add_argument('-n', '--nofilterBam', action='store_true', help='If set, the filtering of BAM files will not be perfomed. Default: Process bam files.' )
    parser.add_argument('-s','--singleEnd', action='store_true', help='Single end read mappings. Default:Paired-end')
    args = parser.parse_args()

    for file in args.bam_files:
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
                bitflagList = filterSingleEndBamAlignments(file,args.o,args.m)
                logging.info("Extracting bitFlag from BAM and adding it to score columns in BED file..")
                replaceBedScoreToBitFlag(outBed,bitflagList)
            else:
                bitflagList = filterBamAlignments(file,args.o,args.m)
                logging.info("Extracting bitFlag from BAM and adding it to score columns in BED file..")
                replaceBedScoreToBitFlag(outBed,bitflagList)

        else:
            logging.info("No filtering of BAM files will be done.")
            bitflagList = bamToBedFromBam(file,args.o,args.m)
            logging.info("Extracting bitFlag from BAM and adding it to score columns in BED file..")
            replaceBedScoreToBitFlag(outBed,bitflagList)



if __name__ == "__main__":
    main()