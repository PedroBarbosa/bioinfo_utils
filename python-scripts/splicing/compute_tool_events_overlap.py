import argparse
from collections import Counter, defaultdict
from matplotlib_venn import venn2, venn3
import matplotlib.pyplot as plt
import upsetplot


def process_vasttools(vastfile, geneid=False):
    try:
        with open(vastfile, 'r') as infile:
            infile.readline()
            genes, events = [], []
            for line in infile:
                if line.split('\t')[0] != "":
                    if line.split('\t')[1] == "":
                        genes.append(line.split('\t')[0])
                    elif geneid:
                        genes.append(line.split('\t')[2])
                    else:
                        genes.append(line.split('\t')[1])

                    events.append(line.split('\t')[7])

            vasttools_counter = Counter(genes)
            vastools_event_types = Counter(events)
        infile.close()
        return vasttools_counter, vastools_event_types

    except TypeError:
        return None, None


def process_majiq(voilafile, geneid=False):
    try:
        with open(voilafile, 'r') as infile:
            infile.readline()
            genes, events = [], []
            for line in infile:
                if line.split('\t')[0] != "":
                    genes.append(line.split('\t')[1]) if not geneid else genes.append(line.split('\t')[2])
                    events.append(line.split('\t')[4])

            majiq_counter = Counter(genes)
            majiq_event_types = Counter(events)
        infile.close()
        return majiq_counter, majiq_event_types

    except TypeError:
        return None, None


def process_rmats(rmatsfile, geneid=False):
    try:
        with open(rmatsfile, 'r') as infile:
            infile.readline()
            genes, events = [], []
            for line in infile:
                if line.split('\t')[0] != "":
                    genes.append(line.split('\t')[1]) if not geneid else genes.append(line.split('\t')[2])

                    events.append(line.split('\t')[3])

            rmats_counter = Counter(genes)
            rmats_type_events = Counter(events)
        infile.close()
        return rmats_counter, rmats_type_events

    except TypeError:
        return None, None


def process_psichomics(psichomicsfile, geneid=False):
    """Finish R psichomics pipeline"""
    genes_in_unique_events = []
    unique_events = []
    events = []
    map_event = {'SE': 'ES'}
    try:
        with open(psichomicsfile, 'r') as infile:
            infile.readline()
            for line in infile:
                if line.split('\t')[0] in unique_events:
                    continue
                unique_events.append(line.split('\t')[0])
                genes_in_unique_events.append(line.split('\t')[1])
                ev = line.split('\t')[3]
                events.append(map_event.get(ev, ev))
            psichomics_counter = Counter(genes_in_unique_events)
            events_counter = Counter(events)
        infile.close()
        return psichomics_counter, events_counter

    except TypeError:
        return None, None


def compute_overlaps(d, gene_id = False):
    final_table = defaultdict(list)
    sets = []

    for tool, counter in d.items():
        sets.append((tool, set(counter.keys())))
        for gene in counter:
            final_table[gene].append([tool, counter[gene]])

    with open("tools_overlap.tsv", "w") as outw:

        if gene_id:
            import mygene
            mg = mygene.MyGeneInfo()
            gene_ids = [k for k in final_table.keys()]
            symbols_dict = mg.querymany(qterms=gene_ids, scopes="ensembl.gene", fields=["symbol"],
                                   returnall=True,
                                   as_dataframe=True,
                                   size=1,
                                   species="human")['out'][["symbol"]].to_dict()['symbol']

            outw.write("#gene\tsymbol\ttools_with_event\tnumber_of_events" + "\n")
        else:
            outw.write("#gene\ttools_with_event\tnumber_of_events" + "\n")
        for k, v in final_table.items():
            flat_list = [item for sublist in v for item in sublist]
            tools = [t for t in flat_list if isinstance(t, str)]
            number_of_events = [str(t) for t in flat_list if isinstance(t, int)]
            if gene_id:
                try:
                    outw.write(k + "\t" + symbols_dict[k] + "\t" + ",".join(tools) + "\t" + ",".join(number_of_events) + "\n")
                except TypeError: #e.g. cases where a vast-tools eventID is reported
                    outw.write(k + "\t" + "" + "\t" + ",".join(tools) + "\t" + ",".join(number_of_events) + "\n")
            else:
                outw.write(k + "\t" + ",".join(tools) + "\t" + ",".join(number_of_events) + "\n")

    outw.close()
    plt.figure()
    plt.title("Genes with splicing events")
    plt.tight_layout()
    out = "tools_overlap.pdf"
    labels = [v[0] for v in sets]
    just_sets = [v[1] for v in sets]
    if len(just_sets) == 2:
        venn2(just_sets, tuple(labels))
    elif len(just_sets) > 3:
        data = upsetplot.from_contents(d)
        upset = upsetplot.UpSet(data, subset_size="count", intersection_plot_elements=4, show_counts='%d')
        upset.plot()
    else:
        venn3(just_sets, tuple(labels))

    plt.savefig(out)
    plt.close()


def write_counts(input_dict):
    out = defaultdict(list)
    tool_order = []
    majiq_complex = 0
    d = {k : v for k,v in input_dict.items() if v}
    for i, tool in enumerate(d.keys()):
        if d[tool]:
            tool_order.append(tool)
            for event_type, counts in d[tool].items():
                if ";" in event_type:
                    majiq_complex += counts
                    continue
                elif event_type not in out.keys():
                    out[event_type] = ["-"] * len(d.keys())

                out[event_type][i] = counts

            if tool == "majiq":
                out["Complex"] = ["-"] * len(d.keys())
                out["Complex"][i] = majiq_complex

    with open("events_table.csv", "w") as file:
        file.write('Event_type' + "\t" + "\t".join(tool_order) + "\n")
        for ev, c in out.items():
            file.write(ev + "\t" + "\t".join(map(str, c)) + "\n")


def main():
    parser = argparse.ArgumentParser(description='Script to compute gene-based event overlaps across different '
                                                 'splicing methods.')
    parser.add_argument('-g', '--gene_id', action='store_true', help='Whether overlaps should be done based on' \
                                                                     'ensembl gene IDs, rather than gene names.')
    parser.add_argument('-v', '--vasttools', help='Path to the output file of process_vastools_compare script (which'
                                                  'picks output of vasttools compare utility and filters by any event '
                                                  'with deltaPSI > 0.2 be default')
    parser.add_argument("-m", "--majiq", help='Path to the output file of process_voila_tsv script where all LSVs are'
                                              'presented (if multiple conditions tested, results come processed.')
    parser.add_argument("-r", "--rmats", help='Path to the output file of process_rmats_maser script where all splicing'
                                              ' events are presented (il multiple conditions tested, results come '
                                              'processed')
    parser.add_argument("-p", "--psichomics", help='Path to the output file of psichomics processing pipeline where all'
                                                   'relevant splicing events are presented (il multiple conditions tested,'
                                                   ' results come concatenated.')
    args = parser.parse_args()
    is_given = 0
    for arg in vars(args):
        if getattr(args, arg):
            is_given += 1

    if is_given < 2:
        raise SystemExit("Minimum of two files are required to compute overlaps. {} file(s) were "
                         "given.".format(is_given))

    vast, vast_c = process_vasttools(args.vasttools, args.gene_id)
    rmats, rmats_c = process_rmats(args.rmats, args.gene_id)
    majiq, majiq_c = process_majiq(args.majiq, args.gene_id)
    psichomics, psichomics_c = process_psichomics(args.psichomics, args.gene_id)

    c = {"vast-tools": vast_c, "rMATS": rmats_c, "MAJIQ": majiq_c, "psichomics": psichomics_c}
    write_counts(c)
    d = {"vast-tools": vast, "rMATS": rmats, "MAJIQ": majiq, "psichomics": psichomics}
    filtered = {k: v for k, v in d.items() if v is not None}
    d.clear()
    d.update(filtered)
    compute_overlaps(d, args.gene_id)


if __name__ == "__main__":
    main()
