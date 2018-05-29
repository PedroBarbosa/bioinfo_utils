import argparse
import os.path
from pybedtools import BedTool

def check_name_col(f):
    if not f.name or not f.strand:
        print('Name and strand columns are required to detect what is the first exon of a gene (crucial to get the promotor region)')
        exit(1)

def process_genome(genome):
    d={}
    with open(genome) as f:
        for line in f:
            (key, val) = line.rstrip().split()
            d[key] = val
    return d

def apply_subset(f,l):
    if len(f) <= l*2:
        print("{}\tfeature can't be subset. Too small".format(str(f).rstrip()))
    else:
        f.start = f.start + l
        f.end = f.end - l
        return f

def subset(bed,n,out):
    bedobj = BedTool(bed)
    bedobj.each(apply_subset, n).saveas(out)

def promoter(bed,n,genome,out):
    gnm_d=process_genome(genome)
    bedobj = BedTool(bed)
    bedobj.each(check_name_col)
    newgene=""
    for f in bedobj:
        if f.name != newgene:
            if f.strand == "+":
                f.start=f.start - n
                if f.start < 0:
                    f.start = 0
                BedTool(f).saveas(out)
                newgene=f.name
            elif f.strand == "-":
                a=1

        elif f.strand == "-":
            a=1


                #BedTool('{} {} {}'.format(vcfrecord.CHROM, vcfrecord.POS - 1, vcfrecord.POS), from_string=True)
            f.slop(l=n,s=True,g=genome).saveas(out)

        else:
            f.saveas(out)
        featname=f.name


def main():
    parser = argparse.ArgumentParser(
        description='Script to play around with bed features depending on specific flags one provide')
    parser.add_argument(dest='bed', help='Path to the bed file')
    parser.add_argument(dest='out', help='Path to the output file')
    parser.add_argument('-n', '--number', type=int, required=True, help='Number of base pairs to play')
    parser.add_argument('-g', '--genome', help='Genome size file required for bedtools slop. (Required if -p is set)')
    m = parser.add_mutually_exclusive_group()
    m.add_argument('-s','--subset', action='store_true', help='Subtract a given number of bp on each feature. E.g: Get deep intronic bed.')
    m.add_argument('-p','--promoter', action='store_true',help='Extend exonic bed file to span likely promotor region. This flag will just look at the first exon of each gene.')
    args = parser.parse_args()

    if args.subset:
        subset(args.bed,args.number,args.out)

    if args.promoter and not args.genome:
        print("Please set '--genome' argument")
        exit(1)
    elif args.promoter:
        if not os.path.isfile(args.genome):
            print("Please set a valid genome file")
            exit(1)
        promoter(args.bed,args.number,args.genome,args.out)


if __name__ == "__main__":
    main()
