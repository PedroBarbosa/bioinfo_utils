import argparse
from collections import defaultdict
import os, re
import pandas as pd
import plotnine as p9
import matplotlib.pyplot as plt


def read_groups(infiles, groups=None):
    """
    Maps voila files to specific majiq run
    provided in the groups argument.

    :param list infiles: List of voila files
    :param str groups: File mapping filename (1st col)
    to group ID (2nd col)
    :return dict: Dict mapping input files and a majiq
    run
    """
    if len(infiles) > 1:
        assert groups is not None, "Multiple files/comparisons were provided. " \
                                   "Please set a groups file accordingly."
    if groups:
        d = defaultdict(list)
        with open(groups, 'r') as f:
            for line in f:
                l = line.rstrip()
                items = l.split("\t")
                key, values = items[0], items[1]
                d[key].append(values)
        f.close()
        if len(d) != len(infiles):
            raise SystemExit("Number of files given is different from "
                             "the number put in the groups file. "
                             "({} vs {})".format(len(infiles), len(d)))
        for f in infiles:
            if not f in d.keys():
                raise SystemExit("File {} is not in the groups configuration".format(f))
        return d

    return


def generate_volcano(data: list,
                     dpsi_threshold: float,
                     prob_threshold: float,
                     outbasename: str):
    """
    Draw volcano plots of differential splicing analysis

    :param list data: List with the data and event types to plot
    :param float dpsi_threshold: dPSI threshold to call an event as significant
    :param float prob_threshold: probability threshold to call an event as significant
    :param str outbasename: Basenome for the output
    """
    print("Generating volcano plots.")

    def _draw(df: pd.DataFrame):

        inc = df[(df.dPSI > dpsi_threshold) & (df.prob > prob_threshold)].shape[0]
        exc = df[(df.dPSI < -dpsi_threshold) & (df.prob > prob_threshold)].shape[0]
        annotations = pd.DataFrame({'xpos': [-0.5, 0.5],
                                   'ypos': [prob_threshold - 0.1, prob_threshold - 0.1],
                                   'annotateText': ["N={}".format(exc), "N={}".format(inc)],
                                   'col': ["goldenrod", "skyblue"],
                                   'size': [30, 30]})

        p1 = (p9.ggplot(df, p9.aes(x='dPSI', y='prob')) +
              p9.geom_point(size=1, colour="grey") +
              p9.geom_point(data=df[(df.dPSI > 0.2) & (df.prob > prob_threshold)], size=1, colour="skyblue") +
              p9.geom_point(data=df[(df.dPSI < -0.2) & (df.prob > prob_threshold)], size=1, colour="goldenrod") +
              p9.scale_colour_manual(["goldenrod", "skyblue"]) +
              p9.geom_text(data=annotations, mapping=p9.aes(x='xpos', y='ypos',
                                                            label='annotateText',
                                                            color='col')) +
              p9.theme_bw() +
              p9.theme(text=p9.element_text(size=12), legend_position="none") +
              p9.geom_vline(xintercept=[-dpsi_threshold, dpsi_threshold], size=0.3) +
              p9.geom_hline(yintercept=prob_threshold, size=0.3) +
              p9.xlab("dPSI change") +
              p9.ylab("Probability (dPSI > {})".format(dpsi_threshold))
              )

        return p1

    def _plot_per_type(data: pd.DataFrame, name: str = "majiq"):
        """
        Dataframe to plot the data
        """
        p9_plot = _draw(data)
        p9_plot.save('{}_{}_all_events_volcano.pdf'.format(outbasename, name), verbose=False)
        plt.close()

        for n, group in data.groupby('event_type'):

            p9_plot = _draw(group)
            p9_plot.save('{}_{}_{}_volcano.pdf'.format(outbasename, name, n), verbose=False)
            plt.close()

    df = pd.DataFrame.from_records(data, columns=['dPSI', 'prob', 'gene_id', 'event_type', 'group'])

    df = df[~df.event_type.str.contains(';')]

    df[['dPSI', 'prob']] = df[['dPSI', 'prob']].apply(pd.to_numeric)
    if df.iloc[0]['group'] != "single_run":
        for name, group in df.groupby('group'):

            _plot_per_type(group, name)

        # print(grouped.name())

    else:
        _plot_per_type(df)


def extract_LSV_types(deltapsi_files, voila_files):
    """
    Extract LSVs type from deltaPSI files.
    (In majiq v2.2 this information was
    removed from the tsv file obtained from
    voila)
    :param list deltapsi_files: List of deltapsi files
    :param list voila:files: List of voila files
    :return dict: Per comparison dictionary with the list
    of splicing event types that each LSV contains
    :return dict: Map of index and filename
    """
    if deltapsi_files is None:
        print("Deltapsi files were not provided. LSV type won't be inferred.")
        return None, None

    assert len(deltapsi_files) == len(voila_files), "Number of deltapsi and voila files must be the same."
    final_lsv_type, file_map = {}, {}
    for i, file in enumerate(deltapsi_files):
        name = os.path.basename(file)
        f = open(file, "r")
        next(f)
        lsvs_dict = {}
        for line in f:
            fields = line.rstrip().split("\t")
            lsv_id = fields[1]
            lsv_type = ';'.join(get_relevant_type(fields[8:11]))
            lsvs_dict[lsv_id] = lsv_type
        final_lsv_type[i] = lsvs_dict
        file_map[i] = name
    return final_lsv_type, file_map


def get_relevant_type(lsvtype):
    """
    Get junction types within a LSV
    :param list lsvtype: List with boolean
    columns for specific splicing junctions
    :return list: Types of junctions present
    in a given LSV
    """

    if isinstance(lsvtype, pd.Series):
        lsvtype = lsvtype.tolist()

    d = {0: "A5SS", 1: "A3SS", 2: "ES"}

    return [d[i] for i, v in enumerate(lsvtype) if v == 'True' or v is True]


def get_new_junctions(junctions_flag: list):
    """
    Returns information about the number of junctions
    within an LSV that pass the significance threshold
    and are known or new
    :param list: junctions_flag: List of binary values
    flagging whether junctions within a LSV
    are known or new (1 known, 0 new)
    :return tuple: Tuple with relevant information
    about the novelty of each significant junciton
    """
    known_junctions = junctions_flag.count("1")
    number_new_junctions = junctions_flag.count("0")
    return (known_junctions, False, number_new_junctions) if number_new_junctions == 0 else (known_junctions, True,
                                                                                             number_new_junctions)


def extract_binary_decisions(lsv_type: str, lsv_type_str: str,
                             dPSIs: str,
                             probabilities_dPSIs: str):
    dPSI_inc_jx, prob_dPSI, inc_idx = 0, 0, 0
    # Pure RI events do not need to know
    # source of target info, because its
    # inclusion on the 2nd group vs 1st group
    # is given in the last jx of the LSV
    if lsv_type_str == "RI":

        if len(dPSIs.split(";")) == 2:
            inc_idx = 1
            dPSI_inc_jx = round(float(dPSIs.split(";")[-1]), 2)
            prob_dPSI = round(float(probabilities_dPSIs.split(";")[-1]), 2)

        else:
            print("Problem in IR event in an LSV.")
            exit(1)

    # binary dPSI orientation of 2nd group vs 1st group
    # is extracted based on source|target tags for ES|ALT3SS|ALT5SS
    #################################
    # source exon - most left junction
    # refers to the inclusion
    ################################
    if lsv_type.split("|")[0] == "s":

        # ALTSS inclusion form is always the shortest transcript
        if len(dPSIs.split(";")) == 2:
            if lsv_type_str in ['ES', 'A5SS']:
                inc_idx = 0
                dPSI_inc_jx = round(float(dPSIs.split(";")[0]), 2)
                prob_dPSI = round(float(probabilities_dPSIs.split(";")[0]), 2)

            elif lsv_type_str == "A3SS":
                inc_idx = 1
                dPSI_inc_jx = round(float(dPSIs.split(";")[1]), 2)
                prob_dPSI = round(float(probabilities_dPSIs.split(";")[1]), 2)

        elif len(dPSIs.split(";")) > 2:
            dPSI_inc_jx, prob_dPSI = 0, 0
            if lsv_type_str in ['ES', 'A5SS']:
                inc_idx = 0
                for idx, x in enumerate(dPSIs.split(";")[:-1]):
                    if abs(float(x)) > abs(dPSI_inc_jx):
                        dPSI_inc_jx = round(float(x), 2)
                        prob_dPSI = round(float(probabilities_dPSIs.split(";")[:-1][idx]), 2)
                        if lsv_type_str == 'ES':
                            inc_idx = idx

            elif lsv_type_str == "A3SS":
                inc_idx = 1
                for idx, x in enumerate(dPSIs.split(";")[1:]):
                    if abs(float(x)) > abs(dPSI_inc_jx):
                        dPSI_inc_jx = round(float(x), 2)
                        prob_dPSI = round(float(probabilities_dPSIs.split(";")[1:][idx]), 2)


    ######################
    # target exon - most left junction
    # refers to the skipping
    ######################
    elif lsv_type.split("|")[0] == "t":

        if lsv_type_str in ['A3SS', 'A5SS']:
            print("SHIT. AltSS when main exon is target(t). {}".format(lsv_type))
            exit(1)

        # binary dPSI orientation of 2nd group vs 1st group
        # is extracted based on source|target tags
        elif len(dPSIs.split(";")) == 2:
            inc_idx = 1
            dPSI_inc_jx = round(float(dPSIs.split(";")[1]), 2)
            prob_dPSI = round(float(probabilities_dPSIs.split(";")[1]), 2)

        elif len(dPSIs.split(";")) > 2:
            dPSI_inc_jx, prob_dPSI = 0, 0
            for idx, x in enumerate(dPSIs.split(";")[1:]):

                if abs(float(x)) > abs(dPSI_inc_jx):
                    dPSI_inc_jx = round(float(x), 2)
                    prob_dPSI = round(float(probabilities_dPSIs.split(";")[1:][idx]), 2)
                    inc_idx = idx

    return dPSI_inc_jx, prob_dPSI, inc_idx


def process_voila_files(files, threshold, probability_threshold, groups, lsv_types, dPSI_file_map):
    """
    Process voila files, possible reading groups set from
    multiple majiq runs as well as capturing LSV splicing
    subtypes by reading deltapsi output files
    :param list files: List of voila files(s). Should not be
    filtered in the voila executable so full gene ranks can be
    obtained and later used for GSEA
    :param float threshold: dPSI threshold to consider a LSV as significant
    :param probability_threshold: probability threshold to consider a
    LSV as significantly different between conditions
    :param dict groups: Dictionary with comparison groups (useful if
    multiple majiq runs were performed. Final output file will analyze
    common significant LSVs between each comparison group)
    :param dict lsv_types: If LSV types were extracted from deltapsi files,
    this info will be written in the output file
    :param dict dPSI_file_map: Just a sanity check dict to make sure
    that corresponding voila and deltapsi files were given in the proper
    order
    :return dict: Significant LSVs dict of that will be outputted
    :return tuple: Genes info for downstream functional enrichment
    :return list: LSV info to draw volcano plots
    """
    final_ranks_to_GSEA, final_all_genes_GO, final_pos_GO_posChange, final_pos_GO_negChange = {}, {}, {}, {}
    final_to_volcano = []
    FINAL_LSVs = defaultdict(list)
    count=0
    for i, f in enumerate(files):
        per_file_ranks, per_file_pos_GO_posChange, \
        per_file_pos_GO_negChange, per_file_all_genes_GO = defaultdict(list), set(), set(), set()
        print("Reading {} file.".format(f))
        if isinstance(lsv_types, dict):
            print("Corresponding deltapsi file from where LSV types were extracted: {}".format(dPSI_file_map[i]))
        group = groups[f][0] if groups else "single_run"

        infile = open(f, 'r')
        for line in infile:
            if not line.startswith("#") and not line.startswith("gene_name"):
                l = line.split("\t")
                gname = l[0]
                gid = l[1].split(".")[0]
                lsv_id = l[2]
                dPSIs = l[3]
                probabilities_dPSIs = l[4]
                lsv_type = l[8]
                njunctions = l[9]
                nexons = l[10]
                is_junction_new = l[11].split(";")
                number_known_junctions, has_new_junction, number_new_junctions = get_new_junctions(is_junction_new)
                chrom = l[12]
                strand = l[13]
                lsv_main_exon_coords = lsv_id.split(":")[-1]
                jx_coord = [i for i in l[14].split(";")]
                coord_exons = [i for i in l[15].split(';') if i != 'nan']
                ri_coords = l[16]

                if isinstance(lsv_types, dict):
                    try:
                        lsv_type_str = lsv_types[i][lsv_id]
                    except KeyError:
                        raise ValueError("Please make sure the order of deltapsi and voila files match the same run")
                    if ri_coords:
                        lsv_type_str += "RI" if lsv_type_str == "" else ";RI"
                else:
                    lsv_type_str = "-"

                # if complex event, extract max dPSI
                if ";" in lsv_type_str:
                    dPSI_inc_jx = str(max([float(x) for x in dPSIs.split(";")], key=abs))
                    prob_dPSI = probabilities_dPSIs.split(";")[dPSIs.split(";").index(dPSI_inc_jx)]
                    lsv_changed_feature = ""

                else:
                    dPSI_inc_jx, prob_dPSI, inc_idx = extract_binary_decisions(lsv_type, lsv_type_str,
                                                                               dPSIs, probabilities_dPSIs)

                    if lsv_type_str == "ES":
                        lsv_changed_feature = chrom + ":" + [x for x in coord_exons if x != lsv_main_exon_coords][
                            inc_idx]

                    elif lsv_type_str in ["A5SS", "A3SS"]:
                        lsv_changed_feature = chrom + ":" + coord_exons[inc_idx]

                    elif lsv_type_str == "RI":
                        lsv_changed_feature = chrom + ":" + ri_coords

                final_to_volcano.append([dPSI_inc_jx, prob_dPSI, gid, lsv_type_str, group])
                per_file_ranks[gid].append((prob_dPSI, dPSI_inc_jx))

                if abs(float(dPSI_inc_jx)) < threshold or float(prob_dPSI) < probability_threshold:
                    per_file_all_genes_GO.add(gid)
                    continue

                _aux = [int(x) for x in re.split(";|-", ';'.join(coord_exons)) if x != 'nan']
                lsv_spanning_coord = chrom + ":" + str(min(_aux)) + "-" + str(max(_aux))
                lsv_main_exon_coords = chrom + ":" + lsv_main_exon_coords
                if float(dPSI_inc_jx) > threshold:
                    per_file_pos_GO_posChange.add(gid)
                elif float(dPSI_inc_jx) < -threshold:
                    per_file_pos_GO_negChange.add(gid)

                _final = [gname, gid, strand, lsv_type_str, str(dPSI_inc_jx), str(prob_dPSI),
                          lsv_main_exon_coords, lsv_changed_feature, lsv_spanning_coord,
                          dPSIs, probabilities_dPSIs, njunctions, nexons,
                          str(number_known_junctions), str(has_new_junction),
                          str(number_new_junctions)]

                if groups is not None:
                    _final.append(group)
                FINAL_LSVs[lsv_id].append(_final)

        final_pos_GO_posChange[group] = per_file_pos_GO_posChange
        final_pos_GO_negChange[group] = per_file_pos_GO_negChange
        final_all_genes_GO[group] = per_file_all_genes_GO
        final_ranks_to_GSEA[group] = per_file_ranks
    return FINAL_LSVs, [final_pos_GO_posChange, final_pos_GO_negChange, final_all_genes_GO,
                        final_ranks_to_GSEA], final_to_volcano


def write_output(final_lsvs, outbasename, enrichment_data, groups):
    """
    Write MAJIQ/Voila output files
    :param dict final_lsvs: Significant LSVs to write as an excel-like file
    :param str outbasename: Basename for the output files
    :param list enrichment_data: Tuple with information to do functional
    enrichment analysis
    :param dict groups: Dictionary with comparison groups
    """

    header = ["#lsv_id", "gene_name", "gene_id", "strand", "lsv_type", "deltaPSI", "probability_deltaPSI",
              "lsv_main_exon_coordinates", "lsv_changed_exon|intron", "lsv_spanning_coordinates",
              "all_delta_PSIs_in_lsv",
              "all_prob_deltaPSIs_in_lsv", "total_numb_junction_in_lsv", "numb_exons_involved", "#known_junctions",
              "is_there_any_new_junction", "#new_junctions"]

    if groups:
        header.append("groups_with_lsv")

    group_overlaps = defaultdict(list)
    strand_map = {"-": "minus", "+": "plus"}

    lsvs_f = open(outbasename + "_majiq_lsvs.tsv", 'w')
    sashimi_f = open(outbasename + "_majiq_toGGsashimi.tsv", 'w')
    bed_f = open(outbasename + "_majiq_lsvs.bed", 'w')
    lsvs_f.write('\t'.join(header) + "\n")

    for lsv_id, data in final_lsvs.items():
        # Per LSV data fields:
        # 0: Gene name
        # 1: Gene ID
        # 2: Strand
        # 3: LSV type
        # 4: max dPSI
        # 5: prob maxDPSI
        # 6: source/target exon coord
        # 7: LSV cassette exon
        # 8: LSV spanning coordinates
        #
        # Iterate over LSV, if found in multiple groups
        for i, v in enumerate(data):
            coord = data[0][8]
            l = re.split(':|-', coord)
            if groups is not None:
                group_overlaps[lsv_id].append(v[-1])
                groups_with_lsv = ';'.join(v[-1] for v in data)

            else:
                groups_with_lsv = outbasename

            if i == 0:
                ev_sash = [v[8], v[0], v[3].replace(";", "_"), strand_map[v[2]], groups_with_lsv, lsv_id]
                ev_bed = [l[0], l[1], l[2], data[0][0] + "_" + data[0][3], groups_with_lsv, data[0][2], lsv_id]
                sashimi_f.write('\t'.join(ev_sash) + '\n')
                bed_f.write('\t'.join(ev_bed) + '\n')
                lsvs_f.write(lsv_id + "\t" + '\t'.join(v) + "\n")

            else:
                lsvs_f.write("\t" + '\t'.join(v) + "\n")

    lsvs_f.close()
    sashimi_f.close()
    bed_f.close()

    if groups is not None:
        lsv_group_overlap_f = open(outbasename + "_majiq_group_maps.tsv", "w")
        for k, v in group_overlaps.items():
            lsv_group_overlap_f.write(k + "\t" + ";".join(v) + "\n")
        lsv_group_overlap_f.close()

    grp = list(enrichment_data[0].keys())
    for g in grp:
        out = outbasename
        if g != "single_run":
            out += "_{}".format(g)

        pos_GO = open(out + "_positive_to_GO.txt", "w")
        pos_GO_posChange = open(out + "_positive_to_GO_pos_dPSIChange.txt", "w")
        pos_GO_negChange = open(out + "_positive_to_GO_neg_dPSIChange.txt", "w")
        neg_GO = open(out + "_negative_to_GO.txt", "w")
        gene_ranks = open(out + "_ranks_to_GSEA.txt", "w")

        pos_go = enrichment_data[0][g].union(enrichment_data[1][g])
        pos_GO.write("\n".join(list(pos_go)))
        pos_GO_posChange.write("\n".join(list(enrichment_data[0][g])))
        pos_GO_negChange.write("\n".join(list(enrichment_data[1][g])))
        neg_GO.write("\n".join(list(enrichment_data[2][g].difference(pos_go))))
        gene_ranks.write("{}\t{}\t{}\n".format("gene_id", "probability_max_dPSI", "max_dPSI"))
        for gene, events in enrichment_data[3][g].items():
            for v in events:
                gene_ranks.write("{}\t{}\t{}\n".format(gene, v[0], v[1]))
        pos_GO.close()
        pos_GO_posChange.close()
        pos_GO_negChange.close()
        neg_GO.close()
        gene_ranks.close()


def main():
    parser = argparse.ArgumentParser(description='Script to retrieve significant LSVs from voila runs')
    parser.add_argument(dest='voila', nargs="+", help='Path to the voila file(s)')
    parser.add_argument("-d", "--deltapsi", metavar="", nargs="+", help='Path to the deltapsi tsv output file(s)'
                                                                        ' so LSV type can be extracted. Order of '
                                                                        'files must match with voila argument.')
    parser.add_argument("-o", "--outbasename", required=True, help='Basename to the output file')
    parser.add_argument("-g", "--groups", metavar="",
                        help='Tab delimited file with groups (2nd col) mapping filenames (1st col)')
    parser.add_argument("-t", "--threshold", metavar="", type=float, default=0.2, help='dPSI threshold. Default:0.2')
    parser.add_argument("-p", "--probability", metavar="", type=float, default=0.9,
                        help='probability threshold. Only LSVs with higher probability than this'
                             'value will be kept. Default:0.9')
    args = parser.parse_args()

    groups = read_groups(args.voila, args.groups)
    lsv_types, dPSI_file_map = extract_LSV_types(args.deltapsi, args.voila)
    final_lsvs, enrichment_data, volcano_data = process_voila_files(args.voila,
                                                                    args.threshold,
                                                                    args.probability,
                                                                    groups,
                                                                    lsv_types,
                                                                    dPSI_file_map)

    generate_volcano(volcano_data, args.threshold, args.probability, args.outbasename)
    write_output(final_lsvs, args.outbasename, enrichment_data, groups)


if __name__ == "__main__":
    main()
