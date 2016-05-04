__author__ = 'pedro'

import argparse
import subprocess
from subprocess import Popen, PIPE, STDOUT
import logging
import sys
import os
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')



def filterBamAlignments(bamfile,outDir,mapper):
    logging.info("Filtering %s file using the specific %s flags for unique mapped reads." % (os.path.basename(bamfile),mapper))

    if mapper == 'bwa_aln':
        nUMR0 = subprocess.Popen(["samtools","view","-h","-F260",bamfile], stdout=subprocess.PIPE)
        #nUMR0.wait()
        nUMR1 = subprocess.Popen(["grep","-w","X0:i:1"], stdin=nUMR0.stdout, stdout=subprocess.PIPE)
        nUMR0.stdout.close()
        nUMR2 = subprocess.Popen(["grep","-wc","X1:i:0"], stdin=nUMR1.stdout, stdout=subprocess.PIPE)
        logging.info("%i alignments kept for further analysis." % int(nUMR2.stdout.read()))
        UMR = subprocess.Popen(["grep","-w","X1:i:0"], stdin=nUMR1.stdout, stdout=subprocess.PIPE)
        nUMR1.stdout.close()
        return UMR.stdout.read()

    elif mapper == 'bwa_mem':
        q10 = subprocess.Popen(["samtools","view", "-h", "-F260", "-q10", bamfile], stdout=subprocess.PIPE)
        #q10.wait()
        nq10XA = subprocess.Popen(["grep", "-vc", "XA:"], stdin=q10.stdout, stdout=subprocess.PIPE)
        logging.info("%i alignments kept for further analysis." % int(nq10XA.stdout.read()))
        q10XA = subprocess.Popen(["grep", "-v", "XA:"], stdin=q10.stdout, stdout=subprocess.PIPE)
        q10.stdout.close()
        return q10XA.stdout.read()

    elif mapper == 'bowtie2':

        q10 = subprocess.Popen(["samtools", "view", "-h", "-F260", "-q10", bamfile],stdout=subprocess.PIPE)
        q10XS = subprocess.Popen(["grep", "-v", "XS:"], stdin=q10.stdout, stdout=subprocess.PIPE)
        q10XSbam = subprocess.Popen(["samtools", "view", "-Shb", "-"], stdin=q10XS.stdout, stderr=subprocess.PIPE, stdout=subprocess.PIPE)


        with open(os.path.join(outDir,os.path.basename(bamfile).replace('.bam', '.bed')), 'w+') as out_file:
            bam2bed = subprocess.Popen(["bedtools", "bamtobed", "-cigar", "-i", "stdin"], stdin=q10XSbam.stdout, stdout=out_file)

        #count elements in bed
        #return subprocess.check_output(['wc', '-l'], stdin=ps.stdout,universal_newlines=True)
        bam2bed.communicate()

  #      out_file.close()
  #      q10.stdout.close()
  #      q10XS.stdout.close()
  #      q10XSbam.stdout.close()


    elif mapper == 'star':
        nq255 = subprocess.Popen(["samtools","view", "-c", "-q255", bamfile], stdout=subprocess.PIPE)
        logging.info("%i alignments kept for further analysis." % int(nq255.stdout.read()))
        q255 = subprocess.Popen(["samtools","view","-h", "-q255", bamfile], stdout=subprocess.PIPE)
        return  q255.stdout.read()

def filterSingleEndBamAlignments(bamfile, mapper):
        logging.info("Filtering %s file using the specific %s flags for unique mapped reads." % (os.path.basename(bamfile),mapper))

        if mapper == 'bowtie2':
            nq255 = subprocess.Popen(["samtools","view", "-c", "-q255", bamfile], stdout=subprocess.PIPE)
            nq255.wait()
            logging.info("%i alignments kept for further analysis." % int(nq255.stdout.read()))
            q255 = subprocess.Popen(["samtools","view","-h", "-q255", bamfile], stdout=subprocess.PIPE)
            nq255.stdout.close()
            return  q255.stdout.read()

        else:
            logging.error("Processing of the single end reads mapping is only available for bowtie2. Sorry!")
            exit(1)

def bamToBedFromBam(bamfile, outputDir, mapper):
    logging.info("Converting BAM %s to BED.." % os.path.basename(bamfile))
    with open(os.path.join(outputDir,os.path.basename(bamfile).replace('.bam', '.bed')), 'w+') as out_file:
        if mapper == "star": #if RNA-seq mappings
            bam2bed = subprocess.Popen(["bedtools", "bamtobed", "-cigar", "-split", "-i", bamfile], stdout=out_file)
            bam2bed.wait()
        else:
            bam2bed = subprocess.Popen(["bedtools", "bamtobed", "-cigar", "-i", bamfile], stdout=out_file)
            bam2bed.wait()

def bamToBedFromBuffer(buffer,bamfile,outputDir, mapper):
    logging.info("Converting filtered BAM %s to BED.." % os.path.basename(bamfile))

    with open(os.path.join(outputDir,os.path.basename(bamfile).replace('.bam', '.bed')), 'w+') as out_file:
        if mapper == "star": #if RNA-seq mappings
            bam2bed = subprocess.Popen(["bedtools", "bamtobed", "-cigar", "-split", "-i", "stdin"], stdin=buffer.stdout, stdout=out_file)
            bam2bed.wait()
        else:
            bam2bed = subprocess.Popen(["bedtools", "bamtobed", "-cigar", "-i", "stdin"], stdin=buffer,stdout=out_file)
            bam2bed.wait()

def main():

    parser = argparse.ArgumentParser(description='Script to analyse the genome coverage and orientation of the scaffolds. Requires samtools and bedtools to be on the system path.')
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
        if not args.nofilterBam:
            if args.singleEnd:

                popen_obj = filterSingleEndBamAlignments(file,args.m)
            else:
                popen_obj = filterBamAlignments(file,args.o,args.m)

            #bamToBedFromBuffer(popen_obj,file,args.o,args.m)


        else:
            logging.info("No filtering of BAM files will be done.")
            bamToBedFromBam(file,args.o,args.m)


if __name__ == "__main__":
    main()