__author__ = 'pedro'

import argparse
import subprocess
import logging
import sys
import os
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')



def filterBamAlignments(bamfile,mapper):
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
        nq10XS = subprocess.Popen(["grep", "-vc", "XS:"], stdin=q10.stdout, stdout=subprocess.PIPE)
        logging.info("%i alignments kept for further analysis." % int(nq10XS.stdout.read()))
        q10XS = subprocess.Popen(["grep", "-v", "XS:"], stdin=q10.stdout, stdout=subprocess.PIPE)
        q10.stdout.close()
        return  q10XS.stdout.read()

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

def bamToBedFromBam(bamfile):

def main():

    parser = argparse.ArgumentParser(description='Script to analyse the genome coverage and orientation of the scaffolds. Requires samtools and bedtools to be on the system path.')
    parser.add_argument(dest='bam_files', metavar='bamFiles', nargs='+', help='Bam files to process.')
    parser.add_argument('-m',metavar = '--mapper', required = True, choices=['bwa_aln','bwa_mem','star','bowtie2'], help='Mapper used in the alignments. Available choices: [bwa_aln,bwa_mem,star,bowtie2].')
    parser.add_argument('-g', metavar = '--genomeTable', required = True, help='Tab delimited genome file. Ex:chrom_name    size(bp). Required by bedtools.')
    parser.add_argument('-f', '--nofilterBam', action='store_true', help='If set, the filtering of BAM files will not be perfomed. Defailt: Process bam files.' )
    parser.add_argument('-s','--singleEnd', action='store_true', help='Single end read mappings. Default:Paired-end')
    args = parser.parse_args()

    for file in args.bam_files:
        if os.path.splitext(file)[1] != ".bam":
            logging.error("Only files with .bam extension are allowed. Exiting!")
            exit(1)

    for file in args.bam_files:
        if not args.nofilterBam:
            if args.singleEnd:

                buffer = filterSingleEndBamAlignments(file,args.m)
            else:
                buffer = filterBamAlignments(file,args.m)

            bamToBed(buffer)


        else:
            bamToBed(file)


if __name__ == "__main__":
    main()