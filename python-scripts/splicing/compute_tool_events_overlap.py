import argparse
from collections import Counter, defaultdict
from matplotlib_venn import venn2, venn3
from itertools import chain
import matplotlib.pyplot as plt
import upsetplot
import mygene
import pandas as pd


def build_gene_table(file: str):
    """
    Builds gene map table. Ensembl geneID
    version numbers are removed, if they exist

    :param str file: Input file
    :return dict: Dictionary with geneIDs as keys
    and gene names as values
    """
    ensembl_gene_map = {}
    with open(file, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("ENS"):
                _l = line.split("\t")
                if _l[0] in ensembl_gene_map:
                    raise ValueError("Duplicate entry in gene map file: {}".format(_l[0]))

                ensembl_gene_map[_l[0]] = _l[1]

    f.close()
    if not ensembl_gene_map:
        raise ValueError("Are you sure the gene map file has Ensembl gene IDs in the first column ?")

    return ensembl_gene_map


def process_vasttools(vast_file: str, genes_map: dict, use_gene_names: bool, keep_vast: bool):
    """
    Process vast-tools processed file and return event counts
    
    :param str vast_file: Path to the input vast-tools processed
    file (which is the output of process_vastools_compare script)
    :param dict genes_map: Dictionary mapping ensembl geneIDs
    to gene names. If `None` and use_gene_names is `False, `
    gene IDs are retrieved from main input file
    :param bool use_gene_names: Whether gene names should be used to
    compute counts and overlaps
    :param bool keep_vast: Keep events that do not have geneID and gene 
    names.
    
    :return Counter: Number of vast-tools events per gene
    :return Counter: Number of vast-tools events per event type
    :return list: List with geneIDs or names that were not
    matched against the genes map.
    """
    unmatched = []
    try:
        with open(vast_file, 'r') as infile:
            infile.readline()
            genes, events, unmatched = [], [], []
            for line in infile:
                if line.split('\t')[0] != "":
                    event_id = line.split('\t')[0]
                    gene_name = line.split('\t')[1]
                    gene_id = line.split('\t')[2]
                    event_type = line.split('\t')[7]

                    if use_gene_names:
                        if gene_name == "":
                            if keep_vast:
                                genes.append(event_id)
                                events.append(event_type)
                            else:
                                print("Vast event ID {} does not have a match "
                                      "to a gene name. Event will be discarded".format(event_id))
                                
                        elif genes_map is None:
             
                            genes.append(gene_name)
                            events.append(event_type)
                                
                        else:
                            names = list(genes_map.values())
                            if gene_name in names:
                                genes.append(gene_name)
                                events.append(event_type)
                            else:
                                unmatched.append(gene_name)

                    else:

                        if gene_id == "":
                            if keep_vast:
                                genes.append(event_id)
                                events.append(event_type)
                            else:
                                print("Vast event ID {} does not have a match "
                                      "to a gene ID. Event will be discarded".format(event_id))

                        elif genes_map is None:
                            genes.append(gene_id)
                            events.append(event_type)

                        else:
                            if gene_id in genes_map:
                                genes.append(gene_id)
                                events.append(event_type)
                            else:
                                unmatched.append(gene_id)

        infile.close()
        vasttools_counter = Counter(genes)
        vastools_event_types = Counter(events)
        return vasttools_counter, vastools_event_types, unmatched

    except TypeError:
        return None, None, None


def process_majiq(voila_file, genes_map, use_gene_names: bool):
    """
    Process MAJIQ processed file and return event counts

    :param str voila_file: Path to the input MAJIQ file (which is the
    output of process_voila_tsv script)
    :param dict genes_map: Dictionary mapping ensembl geneIDs
    to gene names. If `None` and use_gene_names is `False, `
    gene IDs are retrieved from main input file
    :param bool use_gene_names: Whether gene names should be used to
    compute counts and overlaps

    :return Counter: Number of MAJIQ events per gene
    :return Counter: Number of MAJIQ events per event type
    :return list: List with geneIDs or names that were not
    matched against the genes map.
    """
    unmatched = []
    try:
        with open(voila_file, 'r') as infile:
            infile.readline()
            genes, events = [], []

            for line in infile:
                if line.split('\t')[0] != "":
                    event_type = line.split('\t')[4]
                    if event_type:
                        if use_gene_names:
                            gene_name = line.split('\t')[1]

                            if genes_map is None:
                                genes.append(gene_name)
                                events.append(event_type)
                            else:
                                names = list(genes_map.values())
                                if gene_name in names:
                                    genes.append(gene_name)
                                    events.append(event_type)
                                else:
                                    unmatched.append(gene_name)

                        else:
                            gene_id = line.split('\t')[2].split(".")[0]
                            if genes_map is None:
                                genes.append(gene_id)
                                events.append(event_type)

                            else:
                                if gene_id in genes_map:
                                    genes.append(gene_id)
                                    events.append(event_type)
                                else:
                                    unmatched.append(gene_id)

        infile.close()
        majiq_counter = Counter(genes)
        majiq_event_types = Counter(events)

        return majiq_counter, majiq_event_types, unmatched

    except TypeError:
        return None, None, None


def process_rmats(rmats_file: str, genes_map: dict, use_gene_names: bool):
    """
    Process rMATS processed file and return event counts

    :param str rmats_file: Path to the input rMATS file (which is the
    output of process_rmats_maser script)
    :param dict genes_map: Dictionary mapping ensembl geneIDs
    to gene names. If `None` and use_gene_names is `False, `
    gene IDs are retrieved from main input file
    :param bool use_gene_names: Whether gene names should be used to
    compute counts and overlaps
    
    :return Counter: Number of rMATS events per gene
    :return Counter: Number of rMATS events per event type
    :return list: List with geneIDs or names that were not 
    matched against the genes map.
    """
    unmatched = []
    try:
        with open(rmats_file, 'r') as infile:
            infile.readline()
            genes, events = [], []
            for line in infile:
                if line.split('\t')[0] != "":
                    event_type = line.split('\t')[3]
                    # Use gene names
                    if use_gene_names:
                        gene_name = line.split('\t')[1]

                        # If gene map, confirms the existance of
                        # the gene name
                        if genes_map is None:
                            genes.append(gene_name)
                            events.append(event_type)

                        else:
                            names = list(genes_map.values())
                            if gene_name in names:
                                genes.append(gene_name)
                                events.append(event_type)
                            else:
                                unmatched.append(gene_name)

                    # Use gene IDs
                    else:
                        gene_id = line.split('\t')[2].split(".")[0]

                        if genes_map is None:
                            genes.append(gene_id)
                            events.append(event_type)

                        else:
                            if gene_id in genes_map:
                                genes.append(gene_id)
                                events.append(event_type)
                            else:
                                unmatched.append(gene_id)

        infile.close()
        rmats_counter = Counter(genes)
        rmats_type_events = Counter(events)
        return rmats_counter, rmats_type_events, unmatched

    except TypeError:
        return None, None, None


# def process_psichomics(psichomicsfile, ensembl_genes_map, use_gene_id=False, gene_id_from_tool=False):
#     unmatched = []
#     genes = []
#     event_ids = []
#     events = []
#     map_event = {'SE': 'ES'}
#     try:
#         with open(psichomicsfile, 'r') as infile:
#             infile.readline()
#             for line in infile:
#                 if line.split('\t')[0] in event_ids:
#                     continue
#                 else:
#                     event_ids.append(line.split('\t')[0])
# 
#                 if use_gene_id and not gene_id_from_tool:
#                     gene_id = ensembl_genes_map[ensembl_genes_map['gene_name'] == line.split('\t')[1]][
#                         ["gene_name", "gene_id"]]. \
#                         drop_duplicates(keep="first")["gene_id"].to_string(index=False)
# 
#                     if "ENSG" in gene_id:
#                         genes.append(line.split('\t')[1])
#                         ev = line.split('\t')[4]
#                         events.append(map_event.get(ev, ev))
#                     else:
#                         unmatched.append(line.split('\t')[1])
# 
#                 elif gene_id_from_tool:
#                     genes.append(line.split('\t')[2])
#                     ev = line.split('\t')[4]
#                     events.append(map_event.get(ev, ev))
#                 else:
#                     genes.append(line.split('\t')[1])
#                     ev = line.split('\t')[4]
#                     events.append(map_event.get(ev, ev))
# 
#             psichomics_counter = Counter(genes)
#             events_counter = Counter(events)
# 
#         infile.close()
#         return psichomics_counter, events_counter, unmatched
# 
#     except TypeError:
#         return None, None, None


def compute_overlaps(d: dict, unmatched: list, genes_map: dict,
                     use_gene_names: bool, outbasename: str):
    """
    Compute tools overlap

    :param dict d: Dict with event counts per gene, per tool
    :param list unmatched: List of unmatched gene IDs/names
    across any tool
    :param dict genes_map: Ensembl gene map
    :param bool use_gene_names: Whether event counts were performed
    based on gene names
    :param str species: Species to consider in the experiment
    :para str outbasename
    """
    final_table = defaultdict(list)
    sets = []

    if use_gene_names:
        outbasename += "_tools_overlap_gene_names" if genes_map is None else "_tools_overlap_valid_gene_names"
    else:
        outbasename += "_tools_overlap_gene_ids" if genes_map is None else "_tools_overlap_valid_gene_ids"

    for tool, counter in d.items():
        sets.append((tool, set(counter.keys())))
        for gene in counter:
            final_table[gene].append([tool, counter[gene]])

    outw = open("{}.tsv".format(outbasename), "w")

    gene_final_map = {}
    if use_gene_names:
        gene_names = [k for k in final_table.keys()]

        if genes_map is not None:

            for gene in gene_names:
                gene_final_map[gene] = list(genes_map.keys())[list(genes_map.values()).index(gene)]

            if len(unmatched) > 0:
                out_unmatched = open("{}_unmatched_genes.tsv".format(outbasename), "w")
                out_unmatched.write("\n".join(unmatched))
                out_unmatched.close()
            outw.write("#gene_name\tgene_id\ttools_with_event\tnumber_of_events" + "\n")
        else:
            outw.write("#gene_name\ttools_with_event\tnumber_of_events" + "\n")

    else:
        gene_ids = [k for k in final_table.keys()]

        if genes_map is not None:

            for gene_id in gene_ids:
                gene_final_map[gene_id] = genes_map[gene_id]

            if len(unmatched) > 0:
                out_unmatched = open("{}_unmatched_genes.tsv".format(outbasename), "w")
                out_unmatched.write("\n".join(unmatched))
                out_unmatched.close()
            outw.write("#gene_id\tgene_name\ttools_with_event\tnumber_of_events" + "\n")

        else:
            outw.write("#gene_id\ttools_with_event\tnumber_of_events" + "\n")

    for k, v in final_table.items():
        flat_list = [item for sublist in v for item in sublist]
        tools = [t for t in flat_list if isinstance(t, str)]
        number_of_events = [str(t) for t in flat_list if isinstance(t, int)]

        if genes_map is not None:
            outline = [k, gene_final_map[k], ",".join(tools), ",".join(number_of_events)]
        else:
            outline = [k, ",".join(tools), ",".join(number_of_events)]
        outw.write('\t'.join(outline) + "\n")
        # try:
        #     outw.write(k + "\t" + ens_map[k] + "\t" + ",".join(tools) + "\t" + ",".join(number_of_events) + "\n")
        # except TypeError:  # e.g. cases where a vast-tools eventID is reported
        #     outw.write(k + "\t" + "" + "\t" + ",".join(tools) + "\t" + ",".join(number_of_events) + "\n")
        # except KeyError:  # e.g. cases where geneID doesn't exist for the given symbol
        #     outw.write(k + "\t" + "" + "\t" + ",".join(tools) + "\t" + ",".join(number_of_events) + "\n")

    outw.close()
    plt.figure()
    plt.title("Genes with splicing events")
    plt.tight_layout()
    out = "{}.pdf".format(outbasename)
    labels = [v[0] for v in sets]
    just_sets = [v[1] for v in sets]

    # # If ensembl gene map was provided in the first place
    # # include only
    # if gene_final_map is not None:
    #     converted_sets = []
    #
    #     for s in just_sets:
    #         new_s = set()
    #         for gene in s:
    #             try:
    #                 new_s.add(gene_final_map[gene])
    #             except KeyError:
    #                 continue
    #         converted_sets.append(new_s)
    #     just_sets = converted_sets

    if len(just_sets) == 2:
        venn2(just_sets, tuple(labels))
    elif len(just_sets) > 3:

        d = {tool: just_sets[i] for i, tool in enumerate(labels)}
        data = upsetplot.from_contents(d)

        upset = upsetplot.UpSet(data, subset_size="count", intersection_plot_elements=4, show_counts='%d')
        upset.plot()
    else:
        venn3(just_sets, tuple(labels))

    plt.savefig(out)
    plt.close()


def write_event_type_counts(input_dict, gene_map: dict, use_gene_names: bool, outbasename: str):
    """
    Write counts to file

    :param dict input_dict: Dict with event counts per tool
    :param dict gene_map: Ensembl gene map
    :param bool use_gene_names: Whether counts were obtained from gene names
    :param str outbasename
    """
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

    if use_gene_names:
        outbasename += "_gene_names" if gene_map is None else "_valid_gene_names"
    else:
        outbasename += "_gene_ids" if gene_map is None else "_valid_gene_ids"

    with open("{}_events_table.csv".format(outbasename), "w") as file:
        file.write('Event_type' + "\t" + "\t".join(tool_order) + "\n")
        for ev, c in out.items():
            file.write(ev + "\t" + "\t".join(map(str, c)) + "\n")


def main():
    parser = argparse.ArgumentParser(description='Script to compute gene-based event overlaps across different '
                                                 'splicing methods.')

    parser.add_argument(dest='outbasename', help="output basename")
    parser.add_argument('-e', '--ensembl_gene_map', metavar="",
                        help='Tab delimited file mapping Ensembl gene IDs (1st col) '
                             'to gene names (2nd col). Use Biomart, for example, to '
                             'get this file. If set, gene IDs (or gene names '
                             'if \"--use_gene_names\") obtained by each tool will be '
                             'checked against this mapping. Only valid IDs will be '
                             'considered')

    parser.add_argument('-g', '--use_gene_names', action='store_true',
                        help='Whether venn plots should be done based '
                             'on gene names, rather than ensembl gene '
                             'IDs, used by default')

    parser.add_argument('-k', '--keep_vast', action='store_true',
                        help='Keep vast-tools events that did not have a'
                             'gene name/ID match. Default: events are '
                             'discarded')

    parser.add_argument('-v', '--vastools', metavar="", help='Path to the output file of process_vastools_compare '
                                                             'script (which picks output of vast-tools compare utility '
                                                             'and filters by any event with deltaPSI > 0.2 by default')

    parser.add_argument("-m", "--majiq", metavar="", help='Path to the output file of process_voila_tsv script '
                                                          'where all LSVs are presented (if multiple conditions '
                                                          'tested, results come processed.')

    parser.add_argument("-r", "--rmats", metavar="", help='Path to the output file of process_rmats_maser script '
                                                          'where all splicing events are presented (if multiple '
                                                          'conditions tested, results come processed')

    # parser.add_argument("-p", "--psichomics", metavar="", help='Path to the output file of psichomics processing '
    #                                                            'pipeline where all relevant splicing events are '
    #                                                            'presented (if multiple conditions tested, results '
    #                                                            'come concatenated.')
    args = parser.parse_args()
    is_given = 0
    for arg in [args.vastools, args.majiq, args.rmats]: # , args.psichomics]:
        if arg:
            is_given += 1

    if is_given < 2:
        raise ValueError("Minimum of two files are required to compute overlaps. {} "
                         "file(s) were given.".format(is_given))

    if args.ensembl_gene_map:
        genes_map = build_gene_table(args.ensembl_gene_map)
    else:
        genes_map = None

    vast, vast_c, vast_u = process_vasttools(args.vastools,
                                             genes_map,
                                             args.use_gene_names,
                                             args.keep_vast)

    rmats, rmats_c, rmats_u = process_rmats(args.rmats,
                                            genes_map,
                                            args.use_gene_names)

    majiq, majiq_c, majiq_u = process_majiq(args.majiq,
                                            genes_map,
                                            args.use_gene_names)

    # psichomics, psichomics_c, psichomics_u = process_psichomics(args.psichomics,
    #                                                             genes_map,
    #                                                             args.use_gene_names)

    u = [x for x in [vast_u, rmats_u, majiq_u] if x]  # , psichomics_u

    unmatched = list(set(list(chain.from_iterable(u))))

    c = {"vast-tools": vast_c, "rMATS": rmats_c, "MAJIQ": majiq_c}  # , "psichomics": psichomics_c}
    write_event_type_counts(c, genes_map, args.use_gene_names, args.outbasename)
    d = {"vast-tools": vast, "rMATS": rmats, "MAJIQ": majiq}  # , "psichomics": psichomics}
    filtered = {k: v for k, v in d.items() if v is not None}
    d.clear()
    d.update(filtered)
    compute_overlaps(d, unmatched, genes_map, args.use_gene_names, args.outbasename)


if __name__ == "__main__":
    main()
