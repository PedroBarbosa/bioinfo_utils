import argparse
from collections import Counter, defaultdict
from matplotlib_venn import venn2, venn3
from itertools import chain
import matplotlib.pyplot as plt
import upsetplot
import pandas as pd


def retrieve_gene_table(species="human"):
    if species == "human":
        try:
            ensembl_genes_map = pd.read_csv("/Users/pbarbosa/MEOCloud/analysis/genome_utilities/hg38/mart_hg38_v33.txt",
                                            error_bad_lines=False,
                                            sep="\t",
                                            header=0,
                                            usecols=[0, 1, 2, 3, 4, 5, 6],
                                            names=["gene_id", "gene_id_version", "transcript_id",
                                                   "transcript_id_version", "chr", "gene_name", "gene_description"])
        except FileNotFoundError:  # if in lobo
            ensembl_genes_map = pd.read_csv(
                "/mnt/nfs/lobo/MCFONSECA-NFS/mcfonseca/shared/genomes/human/hg38/mart_hg38_v33.txt",
                error_bad_lines=False,
                sep="\t",
                header=0,
                usecols=[0, 1, 2, 3, 4, 5, 6],
                names=["gene_id", "gene_id_version", "transcript_id", "transcript_id_version", "chr", "gene_name",
                       "gene_description"])
    elif species == "mouse":
        try:
            ensembl_genes_map = pd.read_csv(
                "/Users/pbarbosa/MEOCloud/analysis/genome_utilities/mm10/mart_mm10_ensemblv100.txt",
                error_bad_lines=False,
                sep="\t",
                header=0,
                usecols=[1, 2, 3, 4, 5, 6, 7, 8],
                names=["gene_id", "gene_id_version", "transcript_id", "transcript_id_version", "strand",
                       "gene_description", "gene_name", "chr"])
        except FileNotFoundError:
            ensembl_genes_map = pd.read_csv(
                "/mnt/nfs/lobo/MCFONSECA-NFS/mcfonseca/shared/genomes/mouse/GRCm38.p6/mart_mm10_ensemblv100.txt",
                error_bad_lines=False,
                sep="\t",
                header=0,
                usecols=[0, 1, 2, 3, 4, 5, 6, 7],
                names=["gene_id", "gene_id_version", "transcript_id", "transcript_id_version", "strand",
                       "gene_description", "gene_name", "chr"])

    return ensembl_genes_map


def process_vasttools(vastfile, ensembl_genes_map, use_gene_id=False, gene_id_from_tool=False, keep_vast=False):
    unmatched = []
    try:
        with open(vastfile, 'r') as infile:
            infile.readline()
            genes, events = [], []
            for line in infile:
                if line.split('\t')[0] != "":
                    if line.split('\t')[1] == "":
                        if keep_vast:
                            genes.append(line.split('\t')[0])
                            events.append(line.split('\t')[7])
                        else:
                            print("Vast event ID {} does not have a match to a gene name.".format(line.split('\t')[0]))

                    elif use_gene_id and not gene_id_from_tool:
                        gene_id = ensembl_genes_map[ensembl_genes_map['gene_name'] == line.split('\t')[1]][
                            ["gene_name", "gene_id"]]. \
                            drop_duplicates(keep="first")["gene_id"].to_string(index=False)

                        if "ENS" in gene_id:
                            genes.append(line.split('\t')[1])
                            events.append(line.split('\t')[7])
                        else:
                            unmatched.append(line.split('\t')[1])

                    elif use_gene_id and gene_id_from_tool:
                        if "ENS" in line.split('\t')[2]:
                            genes.append(line.split('\t')[2])
                            events.append(line.split('\t')[7])

                    else:
                        genes.append(line.split('\t')[1])
                        events.append(line.split('\t')[7])

            vasttools_counter = Counter(genes)
            vastools_event_types = Counter(events)

        infile.close()
        return vasttools_counter, vastools_event_types, unmatched

    except TypeError:
        return None, None, None


def process_majiq(voilafile, ensembl_genes_map, use_gene_id=False, gene_id_from_tool=False):
    unmatched = []
    try:
        with open(voilafile, 'r') as infile:
            infile.readline()
            genes, events = [], []
            for line in infile:
                if line.split('\t')[0] != "":
                    if use_gene_id and not gene_id_from_tool:
                        gene_id = ensembl_genes_map[ensembl_genes_map['gene_name'] == line.split('\t')[1]][
                            ["gene_name", "gene_id"]]. \
                            drop_duplicates(keep="first")["gene_id"].to_string(index=False)

                        if "ENS" in gene_id:
                            genes.append(line.split('\t')[1])
                            events.append(line.split('\t')[4])
                        else:
                            unmatched.append(line.split('\t')[1])

                    elif gene_id_from_tool:
                        genes.append(line.split('\t')[2])
                        events.append(line.split('\t')[4])

                    else:
                        genes.append(line.split('\t')[1])
                        events.append(line.split('\t')[4])

            majiq_counter = Counter(genes)
            majiq_event_types = Counter(events)
        infile.close()
        return majiq_counter, majiq_event_types, unmatched

    except TypeError:
        return None, None, None


def process_rmats(rmatsfile, ensembl_genes_map, use_gene_id=False, gene_id_from_tool=False):
    unmatched = []
    try:
        with open(rmatsfile, 'r') as infile:
            infile.readline()
            genes, events = [], []
            for line in infile:
                if line.split('\t')[0] != "":
                    if use_gene_id and not gene_id_from_tool:
                        gene_id = ensembl_genes_map[ensembl_genes_map['gene_name'] == line.split('\t')[1]][
                            ["gene_name", "gene_id"]]. \
                            drop_duplicates(keep="first")["gene_id"].to_string(index=False)

                        if "ENS" in gene_id:
                            genes.append(line.split('\t')[1])
                            events.append(line.split('\t')[3])
                        else:
                            unmatched.append(line.split('\t')[1])

                    elif gene_id_from_tool:
                        genes.append(line.split('\t')[2])
                        events.append(line.split('\t')[3])
                    else:
                        genes.append(line.split('\t')[1])
                        events.append(line.split('\t')[3])

            rmats_counter = Counter(genes)
            rmats_type_events = Counter(events)
        infile.close()
        return rmats_counter, rmats_type_events, unmatched

    except TypeError:
        return None, None, None


def process_psichomics(psichomicsfile, ensembl_genes_map, use_gene_id=False, gene_id_from_tool=False):
    unmatched = []
    genes = []
    event_ids = []
    events = []
    map_event = {'SE': 'ES'}
    try:
        with open(psichomicsfile, 'r') as infile:
            infile.readline()
            for line in infile:
                if line.split('\t')[0] in event_ids:
                    continue
                else:
                    event_ids.append(line.split('\t')[0])

                if use_gene_id and not gene_id_from_tool:
                    gene_id = ensembl_genes_map[ensembl_genes_map['gene_name'] == line.split('\t')[1]][
                        ["gene_name", "gene_id"]]. \
                        drop_duplicates(keep="first")["gene_id"].to_string(index=False)

                    if "ENSG" in gene_id:
                        genes.append(line.split('\t')[1])
                        ev = line.split('\t')[4]
                        events.append(map_event.get(ev, ev))
                    else:
                        unmatched.append(line.split('\t')[1])

                elif gene_id_from_tool:
                    genes.append(line.split('\t')[2])
                    ev = line.split('\t')[4]
                    events.append(map_event.get(ev, ev))
                else:
                    genes.append(line.split('\t')[1])
                    ev = line.split('\t')[4]
                    events.append(map_event.get(ev, ev))

            psichomics_counter = Counter(genes)
            events_counter = Counter(events)

        infile.close()
        return psichomics_counter, events_counter, unmatched

    except TypeError:
        return None, None, None


def compute_overlaps(d, unmatched, ensembl_genes_map, species, use_gene_id=False, gene_id_from_tool=False):
    final_table = defaultdict(list)
    sets = []

    if use_gene_id and not gene_id_from_tool:
        outbasename = "tools_overlap_gene_id_fetched_from_symbol"
    elif gene_id_from_tool:
        outbasename = "tools_overlap_gene_id_from_tool"
    else:
        outbasename = "tools_overlap_gene_name"

    for tool, counter in d.items():
        sets.append((tool, set(counter.keys())))
        for gene in counter:
            final_table[gene].append([tool, counter[gene]])

    with open("{}.tsv".format(outbasename), "w") as outw:

        if gene_id_from_tool:
            import mygene
            mg = mygene.MyGeneInfo()
            gene_ids = [k for k in final_table.keys()]
            # Get symbols from geneIDs
            ens_map = mg.querymany(qterms=gene_ids, scopes="ensembl.gene", fields=["symbol"],
                                   returnall=True,
                                   as_dataframe=True,
                                   size=1,
                                   species=species)['out'][["symbol"]].to_dict()['symbol']

            ensembl_map = {k for k, v in ens_map.items() if isinstance(v, float)}
            with open("{}_unmatched_genes.tsv".format(outbasename), "w") as nogene:
                nogene.write("\n".join(list(ensembl_map)) + "\n")
            nogene.close()

        else:
            gene_names = [k for k in final_table.keys()]
            # Get gene IDs from symbols
            # ens_map = mg.querymany(qterms=gene_names, scopes="symbol", fields=["ensembl.gene"], returnall=True,
            #             as_dataframe=True, size=1, species="human")['out'][["ensembl.gene"]].to_dict()['ensembl.gene']
            # ensembl_map = {k for k, v in ens_map.items() if isinstance(v, float)}
            # many unknown as well as genes mapping to multiple IDs (due to the use of haplotype regions)

            ens_map = ensembl_genes_map[ensembl_genes_map['gene_name'].isin(gene_names)][["gene_name", "gene_id"]]. \
                drop_duplicates(keep="first").set_index("gene_name").to_dict()['gene_id']

            with open("{}_unmatched_genes.tsv".format(outbasename), "w") as nogene:
                if not use_gene_id:  # if using gene names, unfetched gene names are only checked here
                    unmatched = [v for v in gene_names if v not in list(ens_map.keys())]

                nogene.write("\n".join(list(unmatched)) + "\n")
            nogene.close()

        if gene_id_from_tool:
            outw.write("#gene_id\tgene_name\ttools_with_event\tnumber_of_events" + "\n")
        else:
            outw.write("#gene_name\tgene_id\ttools_with_event\tnumber_of_events" + "\n")

        for k, v in final_table.items():
            flat_list = [item for sublist in v for item in sublist]
            tools = [t for t in flat_list if isinstance(t, str)]
            number_of_events = [str(t) for t in flat_list if isinstance(t, int)]
            try:
                outw.write(k + "\t" + ens_map[k] + "\t" + ",".join(tools) + "\t" + ",".join(number_of_events) + "\n")
            except TypeError:  # e.g. cases where a vast-tools eventID is reported
                outw.write(k + "\t" + "" + "\t" + ",".join(tools) + "\t" + ",".join(number_of_events) + "\n")
            except KeyError:  # e.g. cases where geneID doesn't exist for the given symbol
                outw.write(k + "\t" + "" + "\t" + ",".join(tools) + "\t" + ",".join(number_of_events) + "\n")

    outw.close()
    plt.figure()
    plt.title("Genes with splicing events")
    plt.tight_layout()
    out = "{}.pdf".format(outbasename)
    labels = [v[0] for v in sets]
    just_sets = [v[1] for v in sets]

    if use_gene_id and not gene_id_from_tool:
        converted_sets = []
        for s in just_sets:
            new_s = set()
            for gene in s:
                try:
                    new_s.add(ens_map[gene])
                except KeyError:
                    continue
            converted_sets.append(new_s)
        just_sets = converted_sets

    if len(just_sets) == 2:
        venn2(just_sets, tuple(labels))
    elif len(just_sets) > 3:
        if use_gene_id:
            d = {tool: just_sets[i] for i, tool in enumerate(labels)}
        data = upsetplot.from_contents(d)

        upset = upsetplot.UpSet(data, subset_size="count", intersection_plot_elements=4, show_counts='%d')
        upset.plot()
    else:
        venn3(just_sets, tuple(labels))

    plt.savefig(out)
    plt.close()


def write_counts(input_dict, use_gene_id=False, gene_id_from_tool=False):
    out = defaultdict(list)
    tool_order = []
    majiq_complex = 0
    d = {k: v for k, v in input_dict.items() if v}

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

            if tool == "MAJIQ":
                out["Complex"] = ["-"] * len(d.keys())
                out["Complex"][i] = majiq_complex

    if use_gene_id and not gene_id_from_tool:
        outbasename = "gene_id_fetched_from_symbol"
    elif gene_id_from_tool:
        outbasename = "gene_id_from_tool"
    else:
        outbasename = "gene_name"

    with open("events_table_{}.csv".format(outbasename), "w") as file:
        file.write('Event_type' + "\t" + "\t".join(tool_order) + "\n")
        for ev, c in out.items():
            file.write(ev + "\t" + "\t".join(map(str, c)) + "\n")


def main():
    parser = argparse.ArgumentParser(description='Script to compute gene-based event overlaps across different '
                                                 'splicing methods.')
    parser.add_argument('-g', '--use_gene_id', action='store_true', help='Whether venn plots should be done based on' \
                                                                         'ensembl gene IDs, rather than gene names.')
    parser.add_argument('-t', '--gene_id_from_tool', action='store_true', help='When -g is true, use gene IDs from '
                                                                               'each tool file, rather than fecthing mygene')
    parser.add_argument('-k', '--keep_vast', action='store_true', help='Keep vast-tools events that did not have a'
                                                                       'gene name/ID match. Default: events are '
                                                                       'discarded')
    parser.add_argument("-s", "--species", type=str, default="human", choices=("human", "mouse"),
                        help='Species. Default:human')
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

    if args.gene_id_from_tool and not args.use_gene_id:
        raise SystemExit("Error. '--use_gene_id' approach is required when using '--gene_id_from_tool' argument")

    ensembl_genes_map = retrieve_gene_table(args.species)
    vast, vast_c, vast_u = process_vasttools(args.vasttools, ensembl_genes_map, args.use_gene_id,
                                             args.gene_id_from_tool, args.keep_vast)
    rmats, rmats_c, rmats_u = process_rmats(args.rmats, ensembl_genes_map, args.use_gene_id, args.gene_id_from_tool)
    majiq, majiq_c, majiq_u = process_majiq(args.majiq, ensembl_genes_map, args.use_gene_id, args.gene_id_from_tool)
    psichomics, psichomics_c, psichomics_u = process_psichomics(args.psichomics, ensembl_genes_map, args.use_gene_id,
                                                                args.gene_id_from_tool)

    u = [x for x in [vast_u, rmats_u, majiq_u, psichomics_u] if x]

    unmatched = list(set(list(chain.from_iterable(u))))
    c = {"vast-tools": vast_c, "rMATS": rmats_c, "MAJIQ": majiq_c, "psichomics": psichomics_c}
    write_counts(c, args.use_gene_id, args.gene_id_from_tool)
    d = {"vast-tools": vast, "rMATS": rmats, "MAJIQ": majiq, "psichomics": psichomics}
    filtered = {k: v for k, v in d.items() if v is not None}
    d.clear()
    d.update(filtered)
    compute_overlaps(d, unmatched, ensembl_genes_map, args.species, args.use_gene_id, args.gene_id_from_tool)


if __name__ == "__main__":
    main()
