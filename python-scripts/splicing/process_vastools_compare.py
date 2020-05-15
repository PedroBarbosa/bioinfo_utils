import argparse
import mygene
import re, os
from collections import defaultdict
from process_voila_tsv import read_groups

def read_groups(groups, infiles):

    d, individual_samples = defaultdict(list), defaultdict(list)
    with open(groups, 'r') as f:
        for line in f:
            l = line.rstrip()
            items = l.split("\t")
            key, values = items[0], items[1]
            d[key].append(values)
            if len(items) > 2:
                individual_samples[key].append((items[2], items[3]))
    if len(d) != len(infiles):
        raise SystemExit("Number of files given is different from the number put in the groups file. ({} vs {})".format(
            len(infiles), len(d)))

    for f in infiles:
        if not f in d.keys():
            raise SystemExit("File {} is not in the groups configuration".format(f))

    return d, individual_samples


def process_vastools_files(files, psi_threshold, groups, individual_samples):

    header_output = ["#Event_ID", "Gene", "GeneID", "Strand", "Coordinates", "Spanning_coordinates", "Event_type", "Major_type",
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
                if not "_" in groups[f][0]:
                    print("Please set comparison ID with \"_\" separating each group. E.g. CTRL_vs_CBE")
                    exit(1)
                group1 = groups[f][0].split("_")[0]
                group2 = groups[f][0].split("_")[-1]
                index_group_1 = [i for i, v in enumerate(header) if v == group1][0]
                index_group_2 = [i for i, v in enumerate(header) if v == group2][0]
                if not index_group_1 or not index_group_2:
                    print("At least one of the groups ({}, {}) is not present in the header of the {}"
                          "file.".format(group1, group2, f))
                    print("Header line: {}".format(header))
                    print("Perhaps you are not using the merged analysis. If so, please add the sample names in the"
                          " groups file. 3rd column: samples of group 1 separated with ','. 4th column: samples of"
                          " group2 separated with ','")
                    exit(1)

            else:
                samples_group_1 = individual_samples[f][0][0]
                samples_group_2 = individual_samples[f][0][1]
                indexes_group_1 = [i for i, v in enumerate(header) if v in samples_group_1]
                indexes_group_2 = [i for i, v in enumerate(header) if v in samples_group_2]
                if len(samples_group_1.split(',')) != len(indexes_group_1) or \
                        len(samples_group_2.split(',')) != len(indexes_group_2):
                    print("At least one of the samples provided ({}, {}) is not present in the header of the {}"
                          " file.".format(samples_group_1, samples_group_2, f))
                    print("Header line: {}".format(header))
            i=1
            for line in infile:
                l = line.rstrip().split("\t")
                gname = l[0]
                event_id = l[1]
                coord = l[2]
                length = l[3]
                event_type = l[5]
                coord_fields = [int(item) for sublist in [re.split('-|,|\+|=', elem)
                                for elem in l[4].split(":")[1:]] for item in sublist
                                if item != ""]

                spanning_coord = "{}:{}-{}".format(l[4].split(":")[0], min(coord_fields), max(coord_fields))
                try:
                    simple_coord = [v for v in re.split(':|-', coord) if v != ""]
                    if len(simple_coord) == 3:
                        if min(coord_fields) > int(simple_coord[1]):
                            print("Minimum spanning coord is higher than coordenate of the event itself.")

                        if max(coord_fields) < int(simple_coord[2]):
                            print("Maximum spanning coord is lower than end coordenate of the event itself.")

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
                    sign_events[event_id].append([gname, coord, spanning_coord, event_type, event_type_map[event_type],
                                                  length, psi_first, psi_second, dpsi, groups[f][0]])

    return sign_events, header_output


def write_output(events, header, outbasename, groups):
    mg = mygene.MyGeneInfo()
    d = defaultdict(list)
    strand_map = {-1: "minus", 1: "plus"}
    strand_map_to_bed = {"minus": "-", "plus": "+"}
    genes = [j[0] for i in events.values() for j in i]
    ensembl_map = mg.querymany(qterms=genes, scopes="symbol", fields=["ensembl.gene", "genomic_pos.strand"], returnall=True,
                          as_dataframe = True, size=1, species="human")['out'][["genomic_pos.strand", "ensembl.gene"]].to_dict()

    with open(outbasename + "_vastools.csv", 'w') as out_main_table:
        out_main_table.write('\t'.join(header) + "\n")
        with open(outbasename + "_vastools_toGGsashimi.csv", 'w') as outsash:
            with open(outbasename + "_vastools.bed", "w") as outbed:
                for event_id, data in events.items():
                    for i, v in enumerate(data):
                        d[event_id + ":" + v[0]].append(v[-1])
                        try:
                            strand = strand_map[ensembl_map['genomic_pos.strand'][v[0]]]
                            geneid = ensembl_map['ensembl.gene'][v[0]]

                        except KeyError:
                            strand = ""

                        v.insert(1, geneid)
                        v.insert(2, strand)
                        if len(groups) > 1:
                            outsash.write("{}\t{}\t{}\t{}\t{}\n".format(v[3], v[0] + "_" + event_id, v[6], v[2], v[-1]))

                        else:
                            outsash.write("{}\t{}\t{}\t{}\n".format(v[3], v[0] + "_" + event_id, v[6], v[2]))

                        if i == 0:
                            out_main_table.write(event_id + "\t" + '\t'.join(v) + "\n")
                        else:
                            out_main_table.write("\t" + '\t'.join(v) + "\n")


                    coord = data[0][3]
                    l = re.split(':|-', coord)
                    ev = [l[0], l[1], l[2], data[0][0] + "_" + data[0][6], ';'.join(v[-1] for v in data)]
                    #if data[0][1]:
                    #    ev.append(strand_map_to_bed[data[0][1]])
                    outbed.write('\t'.join(ev) + '\n')

        if len(groups) > 1:
            outfn_sashimi = outbasename + "_vastools_toGGsashimi.csv"
            with open(outbasename + "_events_group_maps.csv", "w") as out_groups_map:
                for k, v in d.items():
                    out_groups_map.write(k + "\t" + ";".join(sorted(v)) + "\n")

            for group in {groups[key][0] for key in groups.keys()}:
                os.system("grep {} {} | cut -f1,2,3,4 > {}.csv".format(group, outfn_sashimi,
                                                                       outfn_sashimi.replace(".csv", "") +
                                                                       "_" + group))
            with open(outfn_sashimi, "r") as sash_to_remove_group:
                with open(outfn_sashimi + "2.csv", 'w') as sash_group_removed:
                    prev = ""
                    for line in sash_to_remove_group:
                        nogroup = '\t'.join(line.rstrip().split()[:-1])
                        if nogroup != prev:
                            sash_group_removed.write(nogroup + "\n")
                        prev = nogroup
                    os.system("mv {} {}".format(outfn_sashimi + "2.csv", outfn_sashimi))

            out_groups_map.close()
            sash_to_remove_group.close()
            sash_group_removed.close()

    outbed.close()
    out_main_table.close()
    outsash.close()
    os.system("sort -V {} | uniq -u > {}_noDup.csv".format(outbasename + "_vastools.bed", "bed"))
    os.system("mv {} {}".format("bed_noDup.csv", outbasename + "_vastools.bed"))


def main():
    parser = argparse.ArgumentParser(description='Script to produce readable vasttools compare significant events '
                                                 'between two(or more)conditions. It picks from already filtered '
                                                 'events.')
    parser.add_argument(dest='vastools_sign', nargs="+", help='Path to the vasttools files (each file represent and '
                                                           'comparison.')
    parser.add_argument("-o", "--outbasename", required=True, help='Basename to the output file.')
    parser.add_argument("-g", "--groups", required=True, help='Tab delimited file with groups mapping filenames.'
                                                              'If there is only one comparison, the existing file'
                                                              'is enough to later make the overlaps. However,'
                                                              'bed and ggsasshimi files will not be produced.'
                                                              '1st column is the filename, 2nd column the name of '
                                                              'the comparison (group). If analysis was merged, only'
                                                              'these two columns are required. If individual samples'
                                                              'were analysid, a 3rd and 4th column is required with'
                                                              'the individual sample names per group splitted by ","')
    parser.add_argument("-t", "--threshold", type=float, default=20, help='dPSI threshold. Default:20')


    args = parser.parse_args()

    groups, individual_samples = read_groups(args.groups, args.vastools_sign)
    sign_events, header = process_vastools_files(args.vastools_sign, args.threshold, groups, individual_samples)
    write_output(sign_events, header, args.outbasename, groups)


if __name__ == "__main__":
    main()
