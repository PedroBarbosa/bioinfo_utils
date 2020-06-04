import argparse
from collections import defaultdict


def process_file(infile):
    unique_events = defaultdict(list)
    with open(infile, 'r') as f:
        for line in f:
            fields = line.rstrip().split()
            key = ':'.join(fields[:-1])
            group = fields[-1].split(";")
            unique_events[key].append(group)
    f.close()

    for k, combinations in unique_events.items():
        unique_group_combination = set([grp for sublist in combinations for grp in sublist])
        print("\t".join(k.split(":")) + "\t" + ";".join(unique_group_combination))
def main():
    parser = argparse.ArgumentParser(
        description='Script to get the proper view of unique overlapped events by fixing splicing analysis with '
                    'multiple group comparisons.')
    parser.add_argument(dest='input_file', help='Path to the intersected file with the group info in the last column')
    #Example of cmd to generate proper input file
    #bedtools intersect -wo -f 1 -a 0_DCM_paper_diff_splicing_events_hg38.bed -b 2_just_ES_uniq.bed |
    #cut -f1,2,3,4,5,6,11 | sort | uniq > 3_DCM_paper_intersect_full_overlap_with_groups.bed

    args = parser.parse_args()

    process_file(args.input_file)

if __name__ == "__main__":
    main()
