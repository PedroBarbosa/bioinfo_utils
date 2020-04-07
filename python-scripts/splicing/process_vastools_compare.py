import argparse
import os
import glob
from collections import defaultdict
from process_voila_tsv import read_groups

def read_groups(groups, infiles):

    d,individual_samples = defaultdict(list),defaultdict(list)
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

    header_output = ["#Event_ID", "Gene", "Coordinates", "Spanning_coordinates", "Length", "psi_first_group",
              "psi_second_group", "dPSI", "groups_with_event"]

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
            for line in infile:
                l = line.rstrip().split("\t")
                gname = l[0]
                event_id = l[1]
                coord = l[2]
                length = l[3]
                event_type = l[5]
                coord_fields = [item for sublist in [elem.split("-") for elem in l[4].split(":")[1:]] for item in sublist]
                spanning_coord = "{}:{}-{}".format(l[4].split(":")[0], min(coord_fields), max(coord_fields))
                print(spanning_coord, coord, event_type, event_id)
                dpsi = l[-1]
                if individual_samples:
                    psi_first = ','.join([l[i] for i in indexes_group_1])
                    psi_second = ','.join([l[i] for i in indexes_group_2])
                else:
                    psi_first = l[index_group_1]
                    psi_second = l[index_group_2]

                if float(dpsi) >= psi_threshold:
                    sign_events[event_id].append([gname, coord, length, psi_first,
                                              psi_second, dpsi, groups[f][0]])

    return sign_events, header_output


def write_output(events, header, outbasename, groups):
    d = defaultdict(list)
    with open(outbasename + "_vastools.csv", 'w') as out:
        out.write('\t'.join(header) + "\n")
        for event_id, data in events.items():
            for i, v in enumerate(data):
                d[event_id + ":" + v[0]].append(v[-1])
                if i == 0:
                    out.write(event_id + "\t" + '\t'.join(v) + "\n")
                else:
                    out.write("\t" + '\t'.join(v) + "\n")

    out.close()
    if groups:
        with open(outbasename + "_events_group_maps.csv", "w") as out:
            for k, v in d.items():
                out.write(k + "\t" + ";".join(sorted(v)) + "\n")

def main():
    parser = argparse.ArgumentParser(description='Script to produce readable vasttools compare significant events '
                                                 'between two(or more)conditions. It picks from already filtered '
                                                 'events.')
    parser.add_argument(dest='vastools_sign', nargs="+", help='Path to the vasttools files (each file represent and '
                                                           'comparison.')
    parser.add_argument("-o", "--outbasename", required=True, help='Basename to the output file.')
    parser.add_argument("-g", "--groups", required=True, help='Tab delimited file with groups mapping filenames.'
                                                              'If there is only one comparison, the existing file'
                                                              'is enough to later make the overlaps.')
    parser.add_argument("-t", "--threshold", type=float, default=0.2, help='dPSI threshold. Default:0.2')


    args = parser.parse_args()
    groups, individual_samples = read_groups(args.groups, args.vastools_sign)
    sign_events, header = process_vastools_files(args.vastools_sign, args.threshold, groups, individual_samples)
    write_output(sign_events, header, args.outbasename, groups)


if __name__ == "__main__":
    main()
