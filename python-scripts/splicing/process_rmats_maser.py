import argparse
import os
import glob
from collections import defaultdict
from process_voila_tsv import read_groups


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
    header = ["#gene_id", "gene_name", "event_type", "deltaPSI", "pvalue", "FDR", "PSI_1", "PSI_2",
              "coordinates"]
    if groups:
        header.append("group")
    sign_events = defaultdict(list)

    for ev in events:
        print("Processing {} events ".format(ev))
        files = glob.glob(os.path.join(dir, "*sign*{}*".format(ev)))
        for f in files:
            with open(f, 'r') as infile:
                next(infile)
                for line in infile:
                    l = line.rstrip().split("\t")
                    geneid = l[1]
                    genesymbol = l[2]
                    pvalue = l[3]
                    fdr = l[4]
                    deltapsi = l[5]
                    psi_1 = l[6]
                    psi_2 = l[7]
                    coordinates = ':'.join(l[8:14])

                    if float(deltapsi) < psi_threshold or float(fdr) > fdr_threshold:

                        continue

                    if groups:
                        sign_events[geneid].append([genesymbol, ev, deltapsi, pvalue, fdr, psi_1, psi_2, coordinates,
                                                    groups[os.path.basename(f)]])
                    else:
                        sign_events[geneid].append([genesymbol, ev, deltapsi, pvalue, fdr, psi_1, psi_2, coordinates])
            infile.close()
    return sign_events, header


def write_output(events, header, outbasename, groups):
    d = defaultdict(list)
    with open(outbasename + "_rmats.csv", 'w') as out:
        out.write('\t'.join(header) + "\n")
        for geneid, data in events.items():
            for i, v in enumerate(data):
                if groups:
                    d[v[0] + ":" + v[1] + ":" + v[-2]].append(v[-1])
                if i == 0:
                    out.write(geneid + "\t" + '\t'.join(v) + "\n")
                else:
                    out.write("\t" + '\t'.join(v) + "\n")

    out.close()
    if groups:
        with open(outbasename + "_events_group_maps.csv", "w") as out:
            for k, v in d.items():
                out.write(k + "\t" + ";".join(sorted(v)) + "\n")


def main():
    parser = argparse.ArgumentParser(description='Script to produce readable rmats significant events between two(or more)'
                                                 'conditions. It picks from already filtered events through maser so'
                                                 'it\'s possible that by changing thresholds you wont see any'
                                                 'differences.')
    parser.add_argument(dest='rmats_sign', help='Path to the rmats files (each file represent a event type and '
                                                           'comparison')
    parser.add_argument("-o", "--outbasename", required=True, help='Basename to the output file')
    parser.add_argument("-g", "--groups", help='Tab delimited file with groups mapping filenames')
    parser.add_argument("-t", "--threshold", type=float, default=0.2, help='dPSI threshold. Default:0.2')
    parser.add_argument("-f", "--fdr", type=float, default=0.05, help='FDR threshold. Only events with lower FDR value'
                                                                      'will be kept. Default: 0.05')

    args = parser.parse_args()

    if os.path.isdir(args.rmats_sign):
        if args.groups:
            groups = read_groups(args.groups, args.rmats_sign)
        else:
            groups = None

        sign_events, header = process_rmats_files(args.rmats_sign, groups, args.threshold, args.fdr)
        write_output(sign_events, header, args.outbasename, groups)

    else:
        raise SystemExit("rmats output directory given ({}) does not exist.".format(args.rmats_dir))


if __name__ == "__main__":
    main()
