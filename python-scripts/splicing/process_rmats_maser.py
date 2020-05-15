import argparse
import os
import glob
from collections import defaultdict
import re
import numpy as np
import math
import matplotlib.pyplot as plt

def read_groups(groups, dir):
    events = ["SE", "MXE", "A3SS", "A5SS", "RI"]
    d = {}
    with open(groups, 'r') as f:
        for line in f:
            l = line.rstrip()
            items = l.split("\t")
            key, values = items[0], items[1]
            d[key] = values

    files_in_dir = []
    for ev in events:
        files_in_dir.append(glob.glob(os.path.join(dir, "*sign*{}*".format(ev))))
    flat_list = [os.path.basename(item) for sublist in files_in_dir for item in sublist]

    if len(d) != len(flat_list):
        raise SystemExit("Number of files given is different from the number put in the groups file. ({} vs {})".format(
            len(flat_list), len(d)))

    for f in flat_list:
        if not f in d.keys():
            raise SystemExit("File {} is not in the groups configuration".format(f))

    return d


def process_rmats_files(dir, groups, psi_threshold, fdr_threshold):
    events = ["SE", "MXE", "A3SS", "A5SS", "RI"]
    header = ["#coordinates_event_id", "gene_name", "gene_id", "event_type", "deltaPSI", "pvalue", "FDR",
              "PSI_1", "PSI_2", "clean_coordinates"]
    event_type_map = {"SE": "ES"}
    if groups:
        header.append("group")
    sign_events, psi_dist = defaultdict(list), defaultdict(list)

    for ev in events:
        print("Processing {} events ".format(ev))
        files = glob.glob(os.path.join(dir, "*sign*{}*".format(ev)))
        for f in files:
            with open(f, 'r') as infile:
                next(infile)
                for line in infile:
                    l = line.rstrip().split("\t")
                    geneid = l[1].split(".")[0]
                    genesymbol = l[2]
                    pvalue = l[3]
                    fdr = l[4]
                    deltapsi = str(-float(l[5])) #rmats by default display the variation of group1 vs group2, as opposed to the other methods
                    psi_1 = l[6]
                    psi_2 = l[7]
                    coordinates = ':'.join(l[8:14])
                    if abs(float(deltapsi)) < psi_threshold or float(fdr) > fdr_threshold:
                        continue

                    clean_coord = get_clean_coordinates(coordinates, ev)
                    ev = event_type_map.get(ev, ev)
                    if groups:
                        sign_events[coordinates].append([genesymbol, geneid, ev, deltapsi, pvalue, fdr, psi_1, psi_2,
                                                         clean_coord, groups[os.path.basename(f)]])
                    else:
                        sign_events[coordinates].append([genesymbol, geneid, ev, deltapsi, pvalue, fdr, psi_1, psi_2,
                                                         clean_coord])
                    psi_dist[ev].append(deltapsi)

            infile.close()
    return sign_events, header, psi_dist


def get_clean_coordinates(coord, event_type, easy_mode=True):
    coord_fields = coord.split(":")
    chrom = coord_fields[0]
    strand = coord_fields[1]

    coord_fields = [item for sublist in [elem.split("-") for elem in coord_fields[2:]] for item in sublist]
    if easy_mode:
        return "{}:{}-{}".format(chrom, min(coord_fields), max(coord_fields))

    else:
        if event_type in ["SE", "RI"]:
            if strand == "+":
                start, end = coord_fields[2], coord_fields[5]
            elif strand == "-":
                start, end = coord_fields[4], coord_fields[3]

        elif event_type in ["A3SS", "A5SS"]:
            if strand == "+":
                start, end = min(coord_fields[0], coord_fields[2]), coord_fields[5]
            elif strand == "-":
                start, end = coord_fields[4], max(coord_fields[1], coord_fields[3])

        elif event_type == "MXE":
            if strand == "+":
                start, end = coord_fields[4], coord_fields[7]
            elif strand == "-":
                start, end = coord_fields[6], coord_fields[5]

        return "{}:{}-{}".format(chrom, start, end)


def plot_distribution(psi_dist, outbasename):
    for event in psi_dist.keys():
        psis = [float(x) for x in psi_dist[event]]
        if len(psis) > 5:
            plt.figure()
            plt.xlim(-1, 1)
            plt.hist(psis) #, bins=bins)
            plt.title('deltaPSI distribution of {} events (N={})'.format(event, len(psis)))
            plt.xlabel('deltaPSI')
            plt.ylabel('Frequency')
            plt.savefig("{}_{}_psiDist.pdf".format(outbasename, event))


def write_output(events, header, outbasename, groups, tobed):
    d = defaultdict(list)
    strand = {"-": "minus", "+": "plus"}
    with open(outbasename + "_rmats.csv", 'w') as out:
        out.write('\t'.join(header) + "\n")
        with open(outbasename + "_rmats_toGGsashimi.csv", 'w') as outsash:
            with open(outbasename + "_rmats.bed", "w") as outbed:
                for coordinates_id, data in events.items():
                    for i, v in enumerate(data):
                        if groups:
                            full_coord = v[0] + ":" + v[2] + ":" + coordinates_id
                            d[full_coord].append(v[-1])
                            outsash.write("{}\t{}\t{}\t{}\t{}\n".format(v[8], v[0], v[2],
                                                                        strand[coordinates_id.split(":")[1]], v[-1]))
                        else:
                            outsash.write(
                                "{}\t{}\t{}\t{}\n".format(v[8], v[0], v[2], strand[coordinates_id.split(":")[1]]))

                        if i == 0:
                            out.write(coordinates_id + "\t" + '\t'.join(v) + "\n")
                        else:
                            out.write("\t" + '\t'.join(v) + "\n")

                    if tobed:
                        coord = data[0][-2] if groups else data[0][-1]
                        l = re.split(':|-', coord)
                        ev = [l[0], l[1], l[2], data[0][0] + "_" + data[0][2]]
                        if groups:
                            ev.append(';'.join(v[-1] for v in data))
                        outbed.write('\t'.join(ev) + '\n')

    out.close()
    outsash.close()
    outbed.close()

    outfn_sashimi = outbasename + "_rmats_toGGsashimi.csv"
    os.system("sort {} | uniq -u > {}_noDup.csv".format(outfn_sashimi, outfn_sashimi))
    os.system("mv {} {}".format(outfn_sashimi + "_noDup.csv", outfn_sashimi))
    os.system("sort -V {} | uniq -u > {}_noDup.csv".format(outbasename + "_rmats.bed", "bed"))
    os.system("mv {} {}".format("bed_noDup.csv", outbasename + "_rmats.bed"))

    if groups:
        with open(outbasename + "_events_group_maps.csv", "w") as out:
            for k, v in d.items():
                out.write(k + "\t" + ";".join(sorted(v)) + "\n")

        for group in set(groups.values()):
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

        sash_to_remove_group.close()
        sash_group_removed.close()
        os.system("mv {} {}".format(outfn_sashimi + "2.csv", outfn_sashimi))


def main():
    parser = argparse.ArgumentParser(
        description='Script to produce readable rmats significant events between two(or more)'
                    'conditions. It picks from already filtered events through maser so'
                    'it\'s possible that by changing thresholds you wont see any'
                    'differences.')
    parser.add_argument(dest='rmats_sign', help='Path to the rmats files (each file represent a event type and '
                                                'comparison')
    parser.add_argument("-o", "--outbasename", required=True, help='Basename to the output file')
    parser.add_argument("-g", "--groups",
                        help='Tab delimited file with groups mapping filenames. File name (1st column),'
                             ' Group name (2nd column)')
    parser.add_argument("-t", "--threshold", type=float, default=0.2, help='dPSI threshold. Default:0.2')
    parser.add_argument("-f", "--fdr", type=float, default=0.05, help='FDR threshold. Only events with lower FDR value'
                                                                      'will be kept. Default: 0.05')
    parser.add_argument("-b", "--bed", action="store_true",
                        help='Also store a bed file with the coordinates spanning the event')
    parser.add_argument("-p", "--plotPSI", action="store_true",
                        help="Whether histograms of PSI distribution should be plotted")

    args = parser.parse_args()

    if os.path.isdir(args.rmats_sign):
        if args.groups:
            groups = read_groups(args.groups, args.rmats_sign)
        else:
            groups = None

        sign_events, header, psi_dist = process_rmats_files(args.rmats_sign, groups, args.threshold, args.fdr)
        if args.plotPSI:
            plot_distribution(psi_dist, args.outbasename)
        write_output(sign_events, header, args.outbasename, groups, args.bed)

    else:
        raise SystemExit("rmats output directory given ({}) does not exist.".format(args.rmats_dir))


if __name__ == "__main__":
    main()
