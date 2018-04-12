import argparse
from _collections import OrderedDict

parser = argparse.ArgumentParser(
    description='Script to increase the number of base-pairs of base-pairs surrounding splice sites. Input file must contain exonic features only')
parser.add_argument(dest='bed', help='Path to the bed file')
parser.add_argument('-n', '--number', type=int, required=True, help='Number of base pairs to slop')
args = parser.parse_args()

with open(args.bed) as inbed:
    dict=OrderedDict()

    for line in inbed:
        print(int(line.split()[1]) - int(before))
        before=line.split()[2]
        fields=line.rstrip().split()
        del fields[3]
        if line.split()[3] not in dict.keys():
            dict[line.split()[3]] = [fields]
        else:
            dict[line.split()[3]].append(fields)

    before = 0
    for k,v in dict.items():
        dict[k][0][2] = int(dict[k][0][2]) + args.number
        dict[k][-1][1] = int(dict[k][-1][1]) - args.number

        for i in range(1,len(v) -1):
            dict[k][i][1] = int(dict[k][i][1]) - args.number
            dict[k][i][2] = int(dict[k][i][2]) + args.number


    for k,v in dict.items():
        for interval in v:
            print(str(interval[0]) + "\t" + str(interval[1]) + "\t" + str(interval[2]) + "\t" + k + "\t" + str(interval[3]) + "\t" + str(interval[4]))

