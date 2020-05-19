import argparse
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3,venn2_unweighted,venn3_unweighted


parser = argparse.ArgumentParser(description='Script to generate venn diagrams given 2 or 3 sets of files.')
parser.add_argument(dest='inputfiles', metavar='inputfiles', nargs='+', help='Each one representing one set. [2 or 3 files].')
parser.add_argument("-o", "--output", required=True, help='Output basename to write the figure.')
args=parser.parse_args()

if len(args.inputfiles) < 2 or len(args.inputfiles) > 3:
    print("Please provide at least 2 sets to compare.")
    exit(1)

list_sets=[]
labels=[]
for file in args.inputfiles:
    s=set()
    l=os.path.basename(file).split(".")[0]
    labels.append(l)
    with open(file, 'r') as infile:
        for line in infile:
            s.add(line.rstrip())
    infile.close()
    list_sets.append(s)

if len(list_sets) == 2:
    plt.figure(figsize=(7,4))
    venn2_unweighted(list_sets,labels)
elif len(list_sets) == 3:
    plt.figure(figsize=(9,5))
    venn3_unweighted(list_sets, labels)
plt.savefig(args.output + ".png")
