import argparse
from collections import defaultdict
import os, re


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
    d = {0: "A5SS", 1: "A3SS", 2: "ES"}
    return [d[i] for i, v in enumerate(lsvtype) if v == 'True']


def get_new_junctions(junctions_flag):
    """
    Returns information about the number of junctions
    within an LSV that pass the significance threshold
    and are known or new
    :param list: junctions_flag: List of binary values
    flagging whether significant junctions within a LSV
    are known or new (1 known, 0 new)
    :return tuple: Tuple with relevant information
    about the novelty of each significant junciton
    """
    known_junctions = junctions_flag.count("1")
    number_new_junctions = junctions_flag.count("0")
    return (known_junctions, False, number_new_junctions) if number_new_junctions == 0 else (known_junctions, True,
                                                                                             number_new_junctions)


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
    """
    final_ranks_to_GSEA, final_pos_GO, final_neg_GO = {}, {}, {}
    FINAL_LSVs = defaultdict(list)

    for i, f in enumerate(files):
        per_file_ranks, per_file_pos_GO, per_file_neg_GO = defaultdict(list), set(), set()
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
                idx_delta_relevant = [i for i, v in enumerate(dPSIs.split(";")) if abs(float(v)) > threshold and
                                      float(probabilities_dPSIs.split(";")[i]) > probability_threshold]

                max_dPSI = str(max([float(x) for x in dPSIs.split(";")], key=abs))
                prob_max = probabilities_dPSIs.split(";")[dPSIs.split(";").index(max_dPSI)]

                per_file_ranks[gid].append((prob_max, max_dPSI))
                if abs(float(max_dPSI)) < threshold or float(prob_max) < probability_threshold:
                    per_file_neg_GO.add(gid)
                    continue

                per_file_pos_GO.add(gid)
                njunctions = l[9]
                nexons = l[10]
                is_junction_new = [l[11].split(";")[i] for i in idx_delta_relevant]
                number_known_junctions, has_new_junction, number_new_junctions = get_new_junctions(is_junction_new)
                chrom = l[12]
                strand = l[13]
                lsv_main_exon_coords = chrom + ":" + lsv_id.split(":")[-1]
                coord_exons = [int(i) for i in re.split(';|-', l[15]) if i != 'nan']
                lsv_spanning_coord = chrom + ":" + str(min(coord_exons)) + "-" + str(max(coord_exons))
                ri_coords = l[16]
                if isinstance(lsv_types, dict):
                    try:
                        lsv_type = lsv_types[i][lsv_id]
                    except KeyError:
                        raise ValueError("Please make sure the order of deltapsi and voila files match the same run")
                    if ri_coords:
                        lsv_type += "RI" if lsv_type == "" else ";RI"
                else:
                    lsv_type = "-"

                _final = [gname, gid, strand, lsv_type, max_dPSI, prob_max,
                          lsv_main_exon_coords, lsv_spanning_coord,
                          dPSIs, probabilities_dPSIs, njunctions, nexons,
                          str(number_known_junctions), str(has_new_junction),
                          str(number_new_junctions)]
                if groups is not None:
                    _final.append(group)
                FINAL_LSVs[lsv_id].append(_final)

        final_pos_GO[group] = per_file_pos_GO
        final_neg_GO[group] = per_file_neg_GO
        final_ranks_to_GSEA[group] = per_file_ranks
    return FINAL_LSVs, [final_pos_GO, final_neg_GO, final_ranks_to_GSEA]


def write_output(final_lsvs, outbasename, enrichment_data, groups):
    """
    Write MAJIQ/Voila output files
    :param dict final_lsvs: Significant LSVs to write as an excel-like file
    :param str outbasename: Basename for the output files
    :param list enrichment_data: Tuple with information to do functional
    enrichment analysis
    :param dict groups: Dictionary with comparison groups
    """

    header = ["#lsv_id", "gene_name", "gene_id", "strand", "lsv_type", "max_deltaPSI_observed", "prob_max_deltaPSI",
              "lsv_main_exon_coordinates", "lsv_spanning_coordinates", "delta_PSIs_observed", "prob_deltaPSIs_observed",
              "total_numb_junction_in_lsv", "numb_exons_involved", "#relevant_known_junctions(>threshold)",
              "is_there_any_new_relevant_junction", "#relevant_new_junctions(>threshold)"]

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
        # 7: LSV spanning coordinates
        # ...
        # Iterate over LSV, if found in multiple groups
        for i, v in enumerate(data):
            coord = data[0][7]
            l = re.split(':|-', coord)
            if groups is not None:
                group_overlaps[lsv_id].append(v[-1])
                groups_with_lsv = ';'.join(v[-1] for v in data)

            else:
                groups_with_lsv = outbasename

            if i == 0:
                ev_sash = [v[7], v[0], v[3].replace(";", "_"), strand_map[v[2]], groups_with_lsv, lsv_id]
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
        neg_GO = open(out + "_negative_to_GO.txt", "w")
        gene_ranks = open(out + "_ranks_to_GSEA.txt", "w")

        pos_GO.write("\n".join(list(enrichment_data[0][g])))
        neg_GO.write("\n".join(list(enrichment_data[1][g])))
        gene_ranks.write("{}\t{}\t{}\n".format("gene_id", "probability_max_dPSI", "max_dPSI"))
        for gene, events in enrichment_data[2][g].items():
            for v in events:
                gene_ranks.write("{}\t{}\t{}\n".format(gene, v[0], v[1]))
        pos_GO.close()
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
                        help='Tab delimited file with groups (2snd col):q mapping filenames (1st col)')
    parser.add_argument("-t", "--threshold", metavar="", type=float, default=0.2, help='dPSI threshold. Default:0.2')
    parser.add_argument("-p", "--probability", metavar="", type=float, default=0.9,
                        help='probability threshold. Only LSVs with higher probability than this'
                             'value will be kept. Default:0')
    args = parser.parse_args()

    groups = read_groups(args.voila, args.groups)
    lsv_types, dPSI_file_map = extract_LSV_types(args.deltapsi, args.voila)
    final_lsvs, enrichment_data = process_voila_files(args.voila,
                                                      args.threshold,
                                                      args.probability,
                                                      groups,
                                                      lsv_types,
                                                      dPSI_file_map)

    write_output(final_lsvs, args.outbasename, enrichment_data, groups)


if __name__ == "__main__":
    main()
