import argparse
from collections import Counter, defaultdict
from matplotlib_venn import venn2, venn3
import matplotlib.pyplot as plt

def process_vasttools(vastfile):

    try:
        with open(vastfile, 'r') as infile:
            infile.readline()
            genes = [line.split('\t')[0] for line in infile]
            vasttools_counter = Counter(genes)
        infile.close()
        return vasttools_counter

    except TypeError:
        return

def process_majiq(voilafile):

    try:
        with open(voilafile, 'r') as infile:
            infile.readline()
            genes = [line.split('\t')[1] for line in infile if line.split('\t')[0] != ""]
            majiq_counter = Counter(genes)
        infile.close()
        return majiq_counter

    except TypeError:
        return


def process_rmats(rmatsfile):

    try:
        with open(rmatsfile, 'r') as infile:
            infile.readline()
            genes = [line.split('\t')[1] for line in infile if line.split('\t')[0] != ""]
            rmats_counter = Counter(genes)
        infile.close()
        return rmats_counter

    except TypeError:
        return


def compute_overlaps(d):
    final_table = defaultdict(list)
    sets=[]
    for tool, counter in d.items():
        sets.append((tool, set(counter.keys())))
        for gene in counter:
            final_table[gene].append([tool, counter[gene]])

    print("#gene\ttools_with_event\tnumber_of_events")
    for k, v in final_table.items():
        flat_list = [item for sublist in v for item in sublist]
        tools = [t for t in flat_list if isinstance(t, str)]
        number_of_events = [str(t) for t in flat_list if isinstance(t, int)]
        print(k + "\t" + ",".join(tools) + "\t" + ",".join(number_of_events))

    plt.figure()
    plt.title("Genes with splicing events")
    plt.tight_layout()
    out = "venn_plot.pdf"
    labels = [v[0] for v in sets]
    just_sets = [v[1] for v in sets]
    if len(just_sets) == 2:
        venn2(just_sets, tuple(labels))
    else:
        venn3(just_sets, tuple(labels))
    plt.savefig(out)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Script to compute gene-based event overlaps across different '
                                                 'splicing methods.')
    parser.add_argument('-v', '--vasttools', help='Path to the output file of vastools inclusion level script (which'
                                                  'picks output of vasttools tidy utility and filters by any event '
                                                  'with deltaPSI > 0.2')
    parser.add_argument("-m", "--majiq", help='Path to the output file of process_voila_tsv script where all LSVs are'
                                              'presented (if multiple conditions tested, results come processed.')
    parser.add_argument("-r", "--rmats", help='Path to the output file of process_rmats_maser script where all splicing'
                                              ' events are presented (il multiple conditions tested, results come '
                                              'processed')
    args = parser.parse_args()
    is_given = 0
    for arg in vars(args):
        if getattr(args, arg):
            is_given += 1

    if is_given < 2:
        raise SystemExit("Minimum of two files are required to compute overlaps. {} file(s) were "
                         "given.".format(is_given))

    vast = process_vasttools(args.vasttools)
    rmats = process_rmats(args.rmats)
    majiq = process_majiq(args.majiq)
    d = {"vasttools" : vast, "rmats" : rmats, "majiq" : majiq}
    filtered = {k: v for k, v in d.items() if v is not None}
    d.clear()
    d.update(filtered)
    compute_overlaps(d)


if __name__ == "__main__":
    main()
