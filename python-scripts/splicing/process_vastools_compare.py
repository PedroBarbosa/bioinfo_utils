import argparse
import pandas as pd
import re, os
import math
from collections import defaultdict
import mygene


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


def read_groups(groups: str, compare_files: str):
    """
    Reads groups from merged or non-merged vastools
    compare procedures.

    :param str groups: Mapping between filenames and
     groups. If `compare_files` refer to a vastools
     compare procedure without merging the sampeles
     before, additional 2 columns must be present
     in the groups file, mapping individual sample
     names on each group.
    :param str compare_files: List of vastools compare
    files

    :return defaultdict: Dict mapping the compare files
     to the group name
    :return defaultdict: Dict mapping the compare files
    to the individual sample names, if analysis if from
    non-merged procedure
    """
    groups_map, individual_samples = defaultdict(list), defaultdict(list)
    with open(groups, 'r') as f:
        for line in f:
            l = line.rstrip()
            items = l.split("\t")
            key, values = items[0], items[1]
            groups_map[key].append(values)
            if len(items) > 2:
                assert len(items) == 4, "If providing non-merged files, 4 column groups file " \
                                        "must be provided."
                individual_samples[key].append((items[2], items[3]))

    if len(groups_map) != len(compare_files):
        raise SystemExit("Number of files given is different from the number "
                         "put in the groups file. ({} vs {})".format(len(infiles), len(groups_map)))

    for f in compare_files:
        if not f in groups_map.keys():
            raise ValueError("File {} is not in the groups configuration".format(f))

    return groups_map, individual_samples


def process_vastools_files(files: list, groups: dict,
                           psi_threshold: int, individual_samples: dict):
    """
    Process vastools files into a readable format

    :param list files: Input file(s) corresponding to a vastools compare run
    :param dict groups: Groups mapping filenames to comparisons
    :param int psi_threshold: Delta PSI threshold to assign an event
    as significant
    :param dict individual_samples: Dict mapping the compare files
    to the individual sample names, if analysis if from
    non-merged procedure

    :return defaultdict: Dictionary with vast-tools events
    that passed dPSI threshold in any of the comparisons
    :return list: List with the output header line
    """
    header_output = ["#Event_ID", "Gene", "GeneID", "Strand", "Coordinates",
                     "Spanning_coordinates", "Event_type", "Major_type",
                     "Length", "psi_first_group", "psi_second_group", "dPSI", "groups_with_event"]

    event_type_map = {"ANN": "ES",
                      "S": "ES",
                      "C1": "ES",
                      "C2": "ES",
                      "C3": "ES",
                      "MIC": "ES",
                      "Alt3": "A3SS",
                      "Alt5": "A5SS",
                      "IR-C": "RI",
                      "IR-S": "RI",
                      "IR": "RI"}

    sign_events = defaultdict(list)
    for f in files:
        print("Reading {} file.".format(f))
        with open(f, 'r') as infile:
            header = infile.readline().rstrip().split("\t")

            if len(individual_samples) == 0:

                if "_" not in groups[f][0]:
                    raise ValueError("Please set comparison ID with \"_\" separating "
                                     "each group. E.g. CTRL_vs_Treatment")

                group1 = groups[f][0].split("_")[0]
                group2 = groups[f][0].split("_")[-1]
                index_group_1 = [i for i, v in enumerate(header) if v == group1][0]
                index_group_2 = [i for i, v in enumerate(header) if v == group2][0]
                if not index_group_1 or not index_group_2:
                    raise ValueError("\nAt least one of the groups ({}, {}) is not "
                                     "present in the header of the {} file. "
                                     "\nHeader line: {}\n Perhaps your analysis refers "
                                     "to a non-merged procedure. If that is the case, "
                                     "please add the sample names in the groups file. "
                                     "3rd column: samples of group 1 separated with ','. "
                                     "4th column: samples of group2 separated with "
                                     "','".format(group1, group2, f, header))

            else:
                samples_group_1 = individual_samples[f][0][0].split(",")
                samples_group_2 = individual_samples[f][0][1].split(",")
                indexes_group_1 = [i for i, v in enumerate(header) if v in samples_group_1]
                indexes_group_2 = [i for i, v in enumerate(header) if v in samples_group_2]

                if len(samples_group_1) != len(indexes_group_1) or len(samples_group_2) != len(indexes_group_2):
                    raise ValueError("\nAt least one of the samples provided ({}, {}) "
                                     "is not present in the header of the {} file.\nHeader "
                                     "line:\n{}".format(samples_group_1, samples_group_2, f, header))

            for line in infile:
                l = line.rstrip().split("\t")
                gname = l[0]
                event_id = l[1]
                target_coord = l[2]
                length = l[3]
                event_type = l[5]
                coord_fields = [int(item) for sublist in [re.split('-|,|\+|=', elem) for elem in l[4].split(":")[1:]]
                                for item in sublist if item != ""]

                # vastools spanning coord seems to go to
                # the boundaries of upstream (end of it)
                # and downstream (beginning of it) exons.
                # so surrounding exons are not included. Check better
                spanning_coord = "{}:{}-{}".format(l[4].split(":")[0], min(coord_fields), max(coord_fields))
                try:
                    simple_coord = [v for v in re.split(':|-', target_coord) if v != ""]

                    if len(simple_coord) == 3:

                        if min(coord_fields) > int(simple_coord[1]):
                            print("Minimum spanning coord is higher than coordenate of the event itself.")
                            print(event_id, gname)
                        if max(coord_fields) < int(simple_coord[2]):
                            print("Maximum spanning coord is lower than end coordenate of the event itself.")
                            print(event_id, gname)

                except ValueError:
                    print("Problem for {} event with the following simple coordinate: {}".format(event_id, coord))
                    continue

                dpsi = l[-1]
                if individual_samples:
                    psi_first = ','.join([l[i] for i in indexes_group_1])
                    psi_second = ','.join([l[i] for i in indexes_group_2])
                else:
                    psi_first = l[index_group_1]
                    psi_second = l[index_group_2]

                if abs(float(dpsi)) >= psi_threshold:
                    sign_events[event_id].append([gname, target_coord, spanning_coord, event_type,
                                                  event_type_map[event_type], length, psi_first,
                                                  psi_second, dpsi, groups[f][0]])

    return sign_events, header_output


def write_output(sign_events: defaultdict, header: list, outbasename: str,
                 species: str, genes_map: dict = None, ):
    """
    Writes output files from processed dictionaries
    
    :param defaultdict sign_events: Events that passed dPSI threshold
    across any comparison
    :param list header: Header of the output file
    :param str outbasename: Output basename
    :param dict groups: Groups
    :param str species: Species of the experiment. It is necessary
    to query mygene services to extract strand of the genes
    :param dict genes_map: Map of ensembl gene IDs to gene names

    :return: 
    """

    groups_events_map = defaultdict(list)
    strand_map_ggshash = {"-": "minus", "+": "plus"}
    strand_map = {-1: "-", 1: "+", math.nan: ""}
    gene_names = [j[0] for i in sign_events.values() for j in i]

    gene_final_map = {}
    absent_gene_names = []
    mg = mygene.MyGeneInfo()

    # If map is provided, gene IDs will be retrieved from here,
    # and strand info will be obtained with myinfo using gene IDs
    # in the queries
    if genes_map is not None:

        for gene in gene_names:
            try:
                gene_final_map[list(genes_map.keys())[list(genes_map.values()).index(gene)]] = gene
            except ValueError:
                absent_gene_names.append(gene)

        strand_info = mg.querymany(qterms=list(gene_final_map.keys()),
                                   scopes="ensemblgene",
                                   fields=["genomic_pos.strand"],
                                   returnall=True, as_dataframe=True,
                                   size=1, species=species)['out']['genomic_pos.strand'].reset_index().\
            drop_duplicates().rename(columns={"query": "gene_id",
                                              "genomic_pos.strand": "strand"})
        strand_info["gene_name"] = strand_info["gene_id"].map(gene_final_map)

    # If map is NOT provided, both gene IDs and strand info
    # will be retrieved with mygene using gene symbols in the
    # queries
    else:
        strand_info = mg.querymany(qterms=gene_names,
                                   scopes="symbol",
                                   fields=["ensembl.gene", "genomic_pos.strand"],
                                   returnall=True, as_dataframe=True,
                                   size=1, species=species)['out'][['ensembl.gene',
                                                                    'genomic_pos.strand']]. \
            drop_duplicates().reset_index().rename(columns={"query": "gene_name",
                                                            "ensembl.gene": "gene_id",
                                                            "genomic_pos.strand": "strand"})

    strand_info = strand_info.replace({'strand': strand_map})

    ###################
    ## Write outputs ##
    ###################
    out_main = open(outbasename + "_vastools.tsv", 'w')
    out_main.write('\t'.join(header) + "\n")

    outsash = open(outbasename + "_vastools_toGGsashimi.csv", 'w')
    bed_spanning = open(outbasename + "_spanning_vastools.bed", "w")
    bed_target = open(outbasename + "_target_vastools.bed", "w")

    for event_id, data in sign_events.items():

        for i, v in enumerate(data):
            groups_events_map[event_id + ":" + v[0]].append(v[-1])

            _df = strand_info[strand_info.gene_name == v[0]]

            # If gene was not found when retrieving strand with myinfo
            if _df.empty:
                # If it was present in the ensembl map provided as a file
                if gene_final_map is not None:
                    try:
                        gene_id = list(gene_final_map.keys())[list(gene_final_map.values()).index(v[0])]
                    except ValueError:
                        gene_id = ""
                else:
                    gene_id = ""
                strand = ""

            else:
                gene_id = "" if _df.gene_id.isnull().all() else _df.gene_id.item()
                strand = _df.strand.item()

            v.insert(1, gene_id)
            v.insert(2, strand)
            groups_with_event = ';'.join(v[-1] for v in data)

            if i == 0:
                out_main.write(event_id + "\t" + '\t'.join(v) + "\n")

                sash = [v[4], v[0] + "_" + event_id, v[6], strand_map_ggshash.get(v[2], v[2]), groups_with_event]
                outsash.write('\t'.join(sash) + '\n')

                _gene = data[0][0] if data[0][0] else event_id
                coord_span = data[0][4]
                _l = re.split(':|-', coord_span)
                ev_bed_span = [_l[0], _l[1], _l[2], _gene + "_" + data[0][6], groups_with_event, v[2]]
                bed_spanning.write('\t'.join(ev_bed_span) + '\n')

                coord_target = data[0][3]
                _l = re.split(':|-', coord_target)
                ev_bed_target = [_l[0], _l[1], _l[2], _gene + "_" + data[0][6], groups_with_event, v[2]]
                bed_target.write('\t'.join(ev_bed_target) + '\n')

            else:
                out_main.write("\t" + '\t'.join(v) + "\n")

    with open(outbasename + "_events_group_maps.csv", "w") as out_groups_map:
        for k, v in groups_events_map.items():
            out_groups_map.write(k + "\t" + ";".join(sorted(v)) + "\n")

    bed_target.close()
    bed_spanning.close()
    out_main.close()
    outsash.close()
    out_groups_map.close()


def main():
    parser = argparse.ArgumentParser(description='Script to produce readable vast-tools compare significant events '
                                                 'between two (or more) conditions.')

    parser.add_argument(dest='vastools_sign', nargs="+", help='Path to the vasttools files (each file represent and '
                                                              'comparison.')
    parser.add_argument('-o', '--outbasename', required=True, help='Basename to the output file.')
    parser.add_argument('-g', '--groups', required=True, help='Tab delimited file with groups mapping filenames. '
                                                              '1st column is the filename, 2nd column is the name '
                                                              'of the comparison (group). If analysis was merged, '
                                                              'only these two columns are required. If individual '
                                                              'samples were analysid, a 3rd and 4th column are '
                                                              'required with the individual sample names per group '
                                                              'splitted by ",". It is required, even if only one '
                                                              'comparison is under study, because of non-merged '
                                                              'analysis')

    parser.add_argument("-t", "--threshold", type=float, default=20, help='dPSI threshold. Default:20')
    parser.add_argument("-s", "--species", metavar="", type=str, default="human", choices=("human", "mouse"),
                        help='Species. Default:human')

    parser.add_argument('-e', '--ensembl_gene_map', metavar="", help='Tab delimited file mapping Ensembl gene '
                                                                     'IDs (1st col) to gene names (2nd col). '
                                                                     'Use Biomart, for example, to get this file. '
                                                                     'This may be useful to ensure gene names/IDs '
                                                                     'are properly dealt. There were issues in '
                                                                     'the past with gene IDs vast-tools return.')

    args = parser.parse_args()

    if args.ensembl_gene_map:
        genes_map = build_gene_table(args.ensembl_gene_map)
    else:
        genes_map = None

    groups, individual_samples = read_groups(args.groups, args.vastools_sign)

    sign_events, header = process_vastools_files(args.vastools_sign,
                                                 groups,
                                                 args.threshold,
                                                 individual_samples)

    write_output(sign_events, header, args.outbasename, args.species, genes_map)


if __name__ == "__main__":
    main()
