import argparse
from collections import defaultdict
import os, re
import xlsxwriter

def read_groups(groups, infiles):

    d = defaultdict(list)
    with open(groups, 'r') as f:
        for line in f:
            l = line.rstrip()
            items = l.split("\t")
            key, values = items[0], items[1]
            d[key].append(values)

    if len(d) != len(infiles):
        raise SystemExit("Number of files given is different from the number put in the groups file. ({} vs {})".format(
            len(infiles), len(d)))

    for f in infiles:
        if not f in d.keys():
            raise SystemExit("File {} is not in the groups configuration".format(f))

    return d


def get_relevant_type(lsvtype):

    d = {0: "A5SS", 1: "A3SS", 2: "ES"}
    return [d[i] for i, v in enumerate(lsvtype) if v == 'True']


def get_new_junctions(junctions_flag):

    known_junctions = junctions_flag.count("1")
    number_new_junctions = junctions_flag.count("0")
    return (known_junctions, False, number_new_junctions) if number_new_junctions == 0 else (known_junctions, True,
                                                                                             number_new_junctions)


def process_voila_files(files, threshold, probability_threshold, groups=None):
    pos_gens_GO, neg_genes_GO = set(), set()
    rank_genes_to_GSEA = defaultdict(list)
    lsv = defaultdict(list)

    header = ["#lsv_id", "gene_name", "gene_id", "strand", "lsv_type", "max_deltaPSI_observed", "prob_max_deltaPSI",
              "lsv_main_exon_coordinates", "lsv_spanning_coordinates", "delta_PSIs_oberved", "prob_deltaPSIs_observed",
              "total_number_junction_in_lsv", "number_exons_involved", "#relevant_known_junctions(>threshold)",
              "is_there_any_new_relevant_junction", "#relevant_new_junctions(>threshold)"]

    if groups:
        header.append("groups_with_lsv")

    for f in files:

        print("Reading {} file.".format(f))
        with open(f, 'r') as infile:
            next(infile)
            for line in infile:
                l = line.split("\t")
                gname = l[0]
                gid = l[1].split(".")[0]
                lsv_id = l[2]
                delta = l[3]
                probability = l[4]
                idx_delta_relevant = [i for i, v in enumerate(delta.split(";")) if abs(float(v)) > threshold and
                                      float(probability.split(";")[i]) > probability_threshold]

                max_delta = str(max([float(x) for x in delta.split(";")], key=abs))
                prob_max = probability.split(";")[delta.split(";").index(max_delta)]

                rank_genes_to_GSEA[gid.split('.')[0]].append((prob_max, max_delta))
                if abs(float(max_delta)) < threshold or float(prob_max) < probability_threshold:
                    neg_genes_GO.add(gid.split('.')[0])
                    continue

                pos_gens_GO.add(gid.split('.')[0])
                lsv_type = ';'.join(get_relevant_type(l[9:12]))
                njunctions = l[12]
                nexons = l[13]
                is_junction_new = [l[14].split(";")[i] for i in idx_delta_relevant]
                number_known_junctions, has_new_junction, number_new_junctions = get_new_junctions(is_junction_new)
                chrom = l[15]
                strand = l[16]
                lsv_main_exon_coords = chrom + ":" + lsv_id.split(":")[-1]
                coord_exons = [int(i) for i in re.split(';|-', l[18]) if i != 'nan']
                lsv_spanning_coord = chrom + ":" + str(min(coord_exons)) + "-" + str(max(coord_exons))
                ri_coords = l[19]
                if ri_coords:
                    lsv_type += "RI" if not lsv_type else ";RI"

                if groups:
                    lsv[lsv_id].append([gname, gid, strand, lsv_type, max_delta, prob_max, lsv_main_exon_coords, lsv_spanning_coord,
                                        delta, probability, njunctions, nexons, str(number_known_junctions),
                                        str(has_new_junction), str(number_new_junctions), groups[f][0]])
                else:
                    lsv[lsv_id].append([gname, gid, strand, lsv_type, max_delta, prob_max, lsv_main_exon_coords, lsv_spanning_coord,
                         delta, probability, njunctions, nexons, str(number_known_junctions),
                         str(has_new_junction), str(number_new_junctions)])

    return lsv, header, [pos_gens_GO, neg_genes_GO, rank_genes_to_GSEA]


def write_output(lsvs, header, outbasename, tobed, enrichment_data, groups=None):
    d = defaultdict(list)
    strand_map = {"-": "minus", "+": "plus"}
    pos_GO = open(outbasename + "_majiq_positive_list_to_GO.txt", "w")
    neg_GO = open(outbasename + "_majiq_negative_list_to_GO.txt", "w")
    pos_GO.write("\n".join(list(enrichment_data[0])))
    neg_GO.write("\n".join(list(enrichment_data[1])))
    gene_ranks = open(outbasename + "_majiq_gene_ranks_to_GSEA.txt", "w")
    gene_ranks.write("{}\t{}\t{}\n".format("gene_id", "probability_max_dPSI", "max_dPSI"))
    for gene, events in enrichment_data[2].items():
        for v in events:
            gene_ranks.write("{}\t{}\t{}\n".format(gene, v[0], v[1]))
    pos_GO.close()
    neg_GO.close()
    gene_ranks.close()
    with open(outbasename + "_majiq_lsvs.csv", 'w') as out:
        out.write('\t'.join(header) + "\n")
        with open(outbasename + "_majiq_toGGsashimi.csv", 'w') as outsash:
            with open(outbasename + "_majiq.bed", "w") as outbed:
                for lsv, data in lsvs.items():
                    for i, v in enumerate(data):
                        if groups:
                            d[lsv].append(v[-1])
                            outsash.write("{}\t{}\t{}\t{}\t{}\n".format(v[7], v[0], v[3].replace(";", "_"), strand_map[v[2]], v[-1]))
                        else:
                            outsash.write("{}\t{}\t{}\t{}\n".format(v[7], v[0], v[3].replace(";", "_"), strand_map[v[2]]))

                        if i == 0:
                            out.write(lsv + "\t" + '\t'.join(v) + "\n")
                        else:
                            out.write("\t" + '\t'.join(v) + "\n")

                    if tobed:
                        coord = data[0][7]
                        l = re.split(':|-', coord)
                        if groups:
                            ev = [l[0], l[1], l[2], data[0][0] + "_" + data[0][3], ';'.join(v[-1] for v in data), data[0][2]]
                        else:
                            ev = [l[0], l[1], l[2], data[0][0] + "_" + data[0][3], outbasename, data[0][2]]
                        outbed.write('\t'.join(ev) + '\n')

    out.close()
    outsash.close()
    outbed.close()

    if groups:
        outfn_sashimi = outbasename + "_majiq_toGGsashimi.csv"
        with open(outbasename + "_majiq_lsvs_group_maps.csv", "w") as out:
            for k,v in d.items():
                out.write(k + "\t" + ";".join(v) + "\n")

        for group in {groups[key][0] for key in groups.keys()}:
            os.system("grep {} {} | cut -f1,2,3,4 > {}.csv".format(group, outfn_sashimi,
                                                                   outfn_sashimi.replace(".csv", "") + "_" + group))

        with open(outfn_sashimi, "r") as sash_to_remove_group:
            with open(outfn_sashimi + "2.csv", 'w') as sash_group_removed:
                prev = ""
                for line in sash_to_remove_group:
                    nogroup = '\t'.join(line.rstrip().split()[:-1])
                    if nogroup != prev:
                        sash_group_removed.write(nogroup + "\n")
                    prev = nogroup
                os.system("mv {} {}".format(outfn_sashimi + "2.csv", outfn_sashimi))
        out.close()
        sash_to_remove_group.close()
        sash_group_removed.close()
    os.system("sort -V {} | uniq -u > {}_noDup.csv".format(outbasename + "_majiq.bed", "bed"))
    os.system("mv {} {}".format("bed_noDup.csv", outbasename + "_majiq.bed"))


def main():
    parser = argparse.ArgumentParser(description='Script to produce readable LSV between two conditions using voila')
    parser.add_argument(dest='voila', nargs="+", help='Path to the voila file(s)')
    parser.add_argument("-o", "--outbasename", required=True, help='Basename to the output file')
    parser.add_argument("-g", "--groups", help='Tab delimited file with groups mapping filenames')
    parser.add_argument("-t", "--threshold", type=float, default=0.2, help='dPSI threshold. Default:0.2')
    parser.add_argument("-p", "--probability", type=float, default=0.9, help='probability threshold. Only LSVs with higher probability than this'
                                                    'value will be kept. Default:0')
    parser.add_argument("-b", "--bed", action="store_true",
                        help='Also store a bed file with the coordinates spanning the event')
    args = parser.parse_args()

    if args.groups:
        groups = read_groups(args.groups, args.voila)
        lsvs, header, enrichment_data = process_voila_files(args.voila, args.threshold, args.probability, groups=args.groups)
        write_output(lsvs, header, args.outbasename, args.bed, enrichment_data, groups=groups)
    else:
        lsvs, header, enrichment_data = process_voila_files(args.voila, args.threshold, args.probability)
        write_output(lsvs, header, args.outbasename, args.bed, enrichment_data)


if __name__ == "__main__":
    main()
