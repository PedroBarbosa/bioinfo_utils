import argparse
from collections import defaultdict
import os

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

    d = {0: "S5SS", 1: "A3SS", 2: "ES"}
    return [d[i] for i, v in enumerate(lsvtype) if v == 'True']


def get_new_junctions(junctions_flag):

    known_junctions = junctions_flag.count("1")
    number_new_junctions = junctions_flag.count("0")
    return (known_junctions, False, number_new_junctions) if number_new_junctions == 0 else (known_junctions, True,
                                                                                             number_new_junctions)


def process_voila_files(files, sample_groups = None, threshold = 0.2):

    lsv = defaultdict(list)
    header = ["#lsv_id", "gene_name", "gene_id", "lsv_type", "max_deltaPSI_observed", "delta_PSIs_oberved",
              "total_number_junction_in_lsv", "number_exons_involved", "#relevant_known_junctions(>threshold)",
              "is_there_any_new_relevant_junction", "#relevant_new_junctions(>threshold)"]
    if sample_groups:
        header.append("groups_with_lsv")

    for f in files:

        print("Reading {} file.".format(f))
        with open(f, 'r') as infile:
            next(infile)
            for line in infile:
                l = line.split("\t")
                gname = l[0]
                gid = l[1]
                lsv_id = l[2]
                delta = l[3]
                idx_delta_relevant = [i for i, v in enumerate(delta.split(";")) if abs(float(v)) > threshold]
                max_delta = str(max([float(x) for x in delta.split(";")], key=abs))
                lsv_type = ';'.join(get_relevant_type(l[9:12]))
                njunctions = l[12]
                nexons = l[13]
                is_junction_new = [l[14].split(";")[i] for i in idx_delta_relevant]
                number_known_junctions, has_new_junction, number_new_junctions = get_new_junctions(is_junction_new)

                if sample_groups:
                    lsv[lsv_id].append([gname, gid, lsv_type, max_delta, delta, njunctions, nexons,
                                        str(number_known_junctions), str(has_new_junction),
                                        str(number_new_junctions), sample_groups[f][0]])
                else:
                    lsv[lsv_id].append([gname, gid, lsv_type, max_delta, delta, njunctions, nexons,
                                        str(number_known_junctions), str(has_new_junction),
                                        str(number_new_junctions)])
    return lsv, header


def write_output(lsvs, header, outbasename, groups=None):
    d = defaultdict(list)
    with open(outbasename + "_lsvs.csv", 'w') as out:
        out.write('\t'.join(header) + "\n")
        for lsv, data in lsvs.items():
            for i, v in enumerate(data):
                if groups:
                    d[lsv].append(v[-1])
                if i == 0:
                    out.write(lsv + "\t" + '\t'.join(v) + "\n")
                else:
                    out.write("\t" + '\t'.join(v) + "\n")
    out.close()
    if groups:
        with open(outbasename + "_lsvs_group_maps.csv", "w") as out:
            for k,v in d.items():
                out.write(k + "\t" + ";".join(v) + "\n")


def main():
    parser = argparse.ArgumentParser(description='Script to produce readable LSV between two conditions using voila')
    parser.add_argument(dest='voila', nargs="+", help='Path to the voila file(s)')
    parser.add_argument("-o", "--outbasename", required=True, help='Basename to the output file')
    parser.add_argument("-g", "--groups", help='Tab delimited file with groups mapping filenames')
    args = parser.parse_args()

    if args.groups:
        groups = read_groups(args.groups, args.voila)
        lsvs, header = process_voila_files(args.voila, sample_groups=groups)
        write_output(lsvs, header, args.outbasename, groups=groups)
    else:
        lsvs, header = process_voila_files(args.voila)
        write_output(lsvs, header, args.outbasename)


if __name__ == "__main__":
    main()
