import argparse
import os.path
from pybedtools import BedTool

def check_name_col(f):
    if not f.name or not f.strand:
        print('Name and strand columns are required to detect what is the first exon of a gene (crucial to get the promotor region)')
        exit(1)

def get_feature_name(fname,delim,index):
    if not delim:
        return fname
    else:
        return fname.split(delim)[index]

def process_genome(genome):
    d={}
    with open(genome) as f:
        for line in f:
            (key, val) = line.rstrip().split()
            d[key] = int(val)
    return d

def add_bp_positive_feature(f,n):
    f.start = f.start - n
    if f.start < 0:
        f.start = 0
    return str(f)

def add_bp_negative_feature(f,n,gnm_d):
    f.end = f.end + n
    if f.end > gnm_d[f.chrom]:
        f.end = gnm_d[f.chrom]
    return str(f)

def apply_subset(f,l, keep_small):
    if len(f) <= l*2:
        if keep_small:
            return f
        else:
            print("{}\tfeature won't be subset. Too small".format(str(f).rstrip()))
    else:
        f.start = f.start + l
        f.end = f.end - l
        return f

def subset(bed,n,keep_small,out):
    print("Applying subset")
    bedobj = BedTool(bed)
    bedobj.each(apply_subset, n, keep_small).saveas(out)

def promoter(bed,n,genome,delimiter,index,out):
    if os.path.exists(out):
        os.remove(out)

    gnm_d=process_genome(genome)
    bedobj = BedTool(bed)
    bedobj.each(check_name_col)
    fnames_list = [f.name.split(delimiter)[index] for f in bedobj]
    previous_g,previous_f, firstLine, obj_idx, nrecords="","",True, 0,len(bedobj)
    with open(out,'w') as outfile:
        for f in bedobj:

            gname = get_feature_name(f.name,delimiter,index)
            if obj_idx == len(fnames_list)-1: #if last line [index checking will fail]
                if gname != previous_g and f.strand == "+":
                    outfile.write(add_bp_positive_feature(f, n))
                elif f.strand == "+":
                    outfile.write(str(f))
                elif f.strand == "-":
                    outfile.write(add_bp_negative_feature(f, n, gnm_d))

            elif gname != previous_g and f.strand == "+":
                outfile.write(add_bp_positive_feature(f,n))

            elif gname != previous_g and obj_idx != nrecords and fnames_list[obj_idx + 1] == gname:  # first feature of multifeature gene in negative strand
                outfile.write(str(f))

            elif gname != previous_g: # if negative gene has only one feature
                outfile.write(add_bp_negative_feature(f,n,gnm_d))

            elif gname == previous_g and f.strand == "+" or obj_idx != nrecords and fnames_list[obj_idx + 1] == gname: # if positive strand or not last feature of negative gene
                outfile.write(str(f))

            elif gname == previous_g and fnames_list[obj_idx + 1] != gname: #if last feature of a negative strand gene
                outfile.write(add_bp_negative_feature(f,n,gnm_d))

            obj_idx += 1
            previous_g = gname
    outfile.close()

def main():
    parser = argparse.ArgumentParser(
        description='Script to play around with bed features depending on specific flags one provide')
    parser.add_argument(dest='bed', help='Path to the bed file')
    parser.add_argument(dest='out', help='Path to the output file')
    parser.add_argument('-n', '--number', type=int, required=True, help='Number of base pairs to play')
    parser.add_argument('-g', '--genome', help='Genome size file required for bedtools slop. (Required if -p is set)')
    parser.add_argument('-d', '--delimiter', help='Delimiter character to split discriminative feature/subfeature names. Default: Use all name columns')
    parser.add_argument("-i", "--positionFeature", type=int, help='Position index (0-based) where gene names are located in the "-c" columns, when splitted by "-d" delimiter')
    parser.add_argument('-k', "--keep_small", action='store_true', help='Flag to keep intervals that can not be subset due to their short length (smaller than n*2). Default: discard.')
    m = parser.add_mutually_exclusive_group()
    m.add_argument('-s','--subset', action='store_true', help='Subtract a given number of bp on each feature. E.g: Get deep intronic bed.')
    m.add_argument('-p','--promoter', action='store_true',help='Extend exonic bed file to span likely promotor region. This flag will just look at the first exon of each gene.')
    args = parser.parse_args()
    if args.delimiter and args.positionFeature is None:
        print("When --delimiter is set, you should also indicate --positionFeature")
        exit(1)
    elif args.delimiter:
        delimiter=args.delimiter
        index=args.positionFeature
    else:
        delimiter=None
        index=None

    if args.subset:
        subset(args.bed,args.number,args.keep_small, args.out)

    if args.promoter and not args.genome:
        print("Please set '--genome' argument")
        exit(1)
    elif args.promoter:
        if not os.path.isfile(args.genome):
            print("Please set a valid genome file")
            exit(1)
        promoter(args.bed,args.number,args.genome, delimiter, index,args.out)

if __name__ == "__main__":
    main()
